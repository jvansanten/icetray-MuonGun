/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/WeightCalculator.h>
#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Track.h>
#include <boost/foreach.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <simclasses/I3MMCTrack.h>
#include <phys-services/I3Calculator.h>
#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace I3MuonGun {

double
WeightCalculator::GetWeight(const I3Particle &axis, const BundleConfiguration &bundlespec) const
{
	// Get the sampling surface for this bundle.
	SamplingSurfaceConstPtr surface = generator_->GetInjectionSurface(axis, bundlespec);
	std::pair<double, double> steps = surface->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return 0.;
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundlespec.size();
	
	double norm = generator_->GetGeneratedEvents(h, coszen, bundlespec);
	double rate = (*flux_)(h, coszen, m)*surface->GetDifferentialArea(coszen)/norm;
	BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundlespec) {
		if (m > 1)
			rate *= (*radius_)(h, coszen, m, pair.first);
		rate *= (*energy_)(h, coszen, m, pair.first, pair.second);
	}
	
	return rate;
}

// Possibly throw-away utility function: "track" muons to a fixed surface using the
// same method as WeightCalculatorModule
std::vector<I3Particle>
GetMuonsAtSurface(I3FramePtr frame, SurfaceConstPtr surface)
{
	std::vector<I3Particle> final_states;
	
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
	if (!mctree)
		log_fatal("I3MCTree missing!");
	if (!mmctracks)
		log_fatal("I3MMCTrackList missing!");
	BOOST_FOREACH(const Track &track, Track::Harvest(*mctree, *mmctracks)) {
		std::pair<double, double> steps =
		    surface->GetIntersection(track.GetPos(), track.GetDir());
		double energy = track.GetEnergy(steps.first);
		if (energy > 0) {
			final_states.push_back(track);
			I3Particle &p = final_states.back();
			p.SetEnergy(energy);
			p.SetPos(track.GetPos(steps.first));
			p.SetTime(track.GetTime(steps.first));
			p.SetLength(track.GetLength()-steps.first);
		}
	}
	
	return final_states;
}

namespace {

namespace ublas = boost::numeric::ublas;
typedef ublas::bounded_vector<float, 3> vector;

vector
make_vector(double x, double y, double z)
{
	vector v(3);
	v[0] = x; v[1] = y; v[2] = z;
	
	return v;
}

vector
make_vector(const I3Direction &dir)
{
	return make_vector(dir.GetX(), dir.GetY(), dir.GetZ());
}

vector
make_vector(const I3Position &dir)
{
	return make_vector(dir.GetX(), dir.GetY(), dir.GetZ());
}

inline vector
operator-(const I3Position &p1, const I3Position &p2)
{
	return make_vector(p1.GetX()-p2.GetX(), p1.GetY()-p2.GetY(), p1.GetZ()-p2.GetZ());
}

inline vector
operator-(const vector &v, const I3Position &p)
{
	return v-make_vector(p);
}

inline double
GetRadius(const I3Particle &axis, const I3Position &pos)
{
	vector r = pos-axis.GetPos();
	double l2 = ublas::inner_prod(make_vector(axis.GetDir()), r);
	
	return sqrt(ublas::inner_prod(r, r) - l2);
}

}

/**
 * @brief Interface between WeightCalculator and IceTray
 *
 * WeightCalculatorModule handles the details of extracting energies and
 * radial offsets of muons from an I3MCTree and MMCTrackList.
 */
class WeightCalculatorModule : public I3Module {
public:
	WeightCalculatorModule(const I3Context &ctx) : I3Module(ctx)
	{
		AddOutBox("OutBox");
		AddParameter("Flux", "", flux_);
		AddParameter("RadialDistribution", "", radius_);
		AddParameter("EnergyDistribution", "", energy_);
		AddParameter("Generator", "", generator_);
	}
	
	void Configure()
	{
		GetParameter("Flux", flux_);
		GetParameter("RadialDistribution", radius_);
		GetParameter("EnergyDistribution", energy_);
		GetParameter("Generator", generator_);
		
		if (!flux_)
			log_fatal("No flux configured!");
		if (!radius_)
			log_fatal("No radial distribution configured!");
		if (!energy_)
			log_fatal("No energy distribution configured!");
		if (!generator_)
			log_fatal("No generator configured!");
	}
	
	void DAQ(I3FramePtr frame)
	{
		// First, harvest the muons in the bundle at their points of injection, storing
		// everything that's necessary to estimate the energy lost up to an arbitrary point
		I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
		I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
		if (!mctree)
			log_fatal("I3MCTree missing!");
		if (!mmctracks)
			log_fatal("I3MMCTrackList missing!");
		const I3MCTree::iterator primary = mctree->begin();
		std::list<Track> tracks = Track::Harvest(*mctree, *mmctracks);
		BundleConfiguration bundlespec;
		BOOST_FOREACH(const Track &track, tracks)
			bundlespec.push_back(std::make_pair(
			    GetRadius(*primary, track.GetPos()), track.GetEnergy()));
		
		// Now, track the bundle to the innermost sampling surface, which may or may
		// not be identical to the injection surface. Since we allow the surface to
		// depend on the bundle configuration, we have to repeat this track/query
		// loop until the surface is stable. 
		
		// We will accept tracks that are reasonably close to the proposed surface
		double tol = 1*I3Units::m;
		SamplingSurfaceConstPtr surface;
		bool stable = false;
		while (!stable) {
			stable = true;
			surface = generator_->GetInjectionSurface(*primary, bundlespec);
			std::list<Track>::iterator track = tracks.begin();
			BundleConfiguration::iterator bspec = bundlespec.begin();
			for ( ; track != tracks.end() && bspec != bundlespec.end(); ) {
				std::pair<double, double> steps =
				    surface->GetIntersection(track->GetPos(), track->GetDir());
				if (!std::isfinite(steps.first)) {
					// This track misses the surface; drop it
					track = tracks.erase(track);
					bspec = bundlespec.erase(bspec);
				} else if (steps.first > tol) {
					// This track will pierce the newly-proposed
					// surface, but is still outside. Track it to the
					// new surface and update its energy and radius.
					double energy = track->GetEnergy(steps.first);
					if (energy > 0) {
						bspec->first = GetRadius(*primary, track->GetPos(steps.first));
						bspec->second = energy;
						track++;
						bspec++;
					} else {
						track = tracks.erase(track);
						bspec = bundlespec.erase(bspec);
					}
					// We had to update the bundle and need to check
					// that the proposed surface is stable in the
					// next iteration.
					stable = false;
				} else {
					// This track is at the proposed surface.
					track++;
					bspec++;
				}
			}
		}
		
		double rate = 0.;
		std::pair<double, double> steps =
		    surface->GetIntersection(primary->GetPos(), primary->GetDir());
		if (bundlespec.size() > 0 && std::isfinite(steps.first)) {
			double h = GetDepth(primary->GetPos().GetZ() + steps.first*primary->GetDir().GetZ());
			double coszen = cos(primary->GetDir().GetZenith());
			unsigned m = bundlespec.size();
			
			double norm = generator_->GetGeneratedEvents(h, coszen, bundlespec);
			rate = (*flux_)(h, coszen, m)*surface->GetDifferentialArea(coszen)/norm;
			BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundlespec) {
				if (m > 1)
					rate *= (*radius_)(h, coszen, m, pair.first);
				rate *= (*energy_)(h, coszen, m, pair.first, pair.second);
			}
		}
		
		frame->Put(GetName(), boost::make_shared<I3Double>(rate));
		PushFrame(frame);
	}
private:
	FluxPtr flux_;
	RadialDistributionPtr radius_;
	EnergyDistributionPtr energy_;
	GenerationProbabilityPtr generator_;
};

}

I3_MODULE(I3MuonGun::WeightCalculatorModule);
