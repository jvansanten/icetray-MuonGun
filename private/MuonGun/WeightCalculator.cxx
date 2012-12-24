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

struct Track {
	Track(const I3MMCTrack &t) : pos(make_vector(t.GetXi(), t.GetYi(), t.GetZi())),
	    dir(make_vector(t.GetI3Particle().GetDir())), energy(t.GetEi()), time(t.GetTi()),
	    length(t.GetI3Particle().GetLength()-ublas::norm_2(pos-t.GetI3Particle().GetPos()))
	{};
	Track(const I3Particle &t) : pos(make_vector(t.GetPos())), dir(make_vector(t.GetDir())),
	    energy(t.GetEnergy()), time(t.GetTime()), length(t.GetLength())
	{};
	double radius(const I3Particle &axis) const
	{
		vector r = pos-axis.GetPos();
		double l2 = ublas::inner_prod(make_vector(axis.GetDir()), r);
		
		return sqrt(ublas::inner_prod(r, r) - l2);
	}
	void advance(double l)
	{
		if (l < 0) {
			return;
		} else if (l >= length) {
			length = 0;
			energy = 0;
		} else {
			pos += l*dir;
			time += l/I3Constants::c;
			length -= l;
			for ( ; (secondaries.first != secondaries.second) &&
			    (secondaries.first->GetTime() <= time); secondaries.first++)
				energy -= secondaries.first->GetEnergy();
		}
	}
	
	I3Position GetPos() const { return I3Position(pos[0], pos[1], pos[2]); }
	I3Direction GetDir() const { return I3Direction(dir[0], dir[1], dir[2]); }
	
	vector pos;
	vector dir;
	double energy;
	double time;
	double length;
	std::pair<I3MCTree::sibling_iterator, I3MCTree::sibling_iterator> secondaries;
};

inline bool
operator!=(const I3MMCTrack &track, const I3Particle &p)
{
	return (track.GetI3Particle().GetMajorID() != p.GetMajorID() ||
	    track.GetI3Particle().GetMinorID() != p.GetMinorID());
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
		I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
		I3MCTree::iterator primary = mctree->begin();
		
		// First, harvest the muons in the bundle at their points of injection,
		// storing iterators over the secondary energy losses
		std::vector<Track> tracks;
		if (primary->GetType() == I3Particle::PPlus || primary->IsNucleus()) {
			// For complete air-shower simulation, harvest radial offsets and
			// energies at the sampling surface from the MMCTrackList
			I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
			if (!mmctracks)
				log_fatal("This appears to be CORSIKA simulation, but I have no MMCTrackList!");
			I3MCTree::iterator p = mctree->begin();
			BOOST_FOREACH(const I3MMCTrack &mmctrack, *mmctracks) {
				while (p != mctree->end() && mmctrack != *p)
					p++;
				if (p != mctree->end()) {
					Track track(mmctrack);
					track.secondaries.first = mctree->begin(p);
					track.secondaries.second = mctree->end(p);
					tracks.push_back(track);
					
					p = mctree->end(p);
				}
			}
		} else {
			// Harvest the direct daughters of the primary.
			for (I3MCTree::iterator p = mctree->begin(primary); p != mctree->end(primary); ) {
				if (p->GetType() != I3Particle::MuMinus) {
					p++;
					continue;
				} else {
					Track track(*p);
					track.secondaries.first = mctree->begin(p);
					track.secondaries.second = mctree->end(p);
					tracks.push_back(track);
					
					p = mctree->end(p);
				}
			}
		}
		
		// Record the state of the bundle at the point of injection
		BundleConfiguration bundlespec;
		BOOST_FOREACH(const Track &track, tracks)
			bundlespec.push_back(std::make_pair(track.radius(*primary), track.energy));
		
		double rate = 0.;
		// Get the sampling surface for this bundle.
		SamplingSurfaceConstPtr surface = generator_->GetInjectionSurface(*primary, bundlespec);
		std::pair<double, double> steps = surface->GetIntersection(primary->GetPos(), primary->GetDir());
		if (std::isfinite(steps.first)) {
			
			// Track bundle to the sampling surface, which may be identical
			// to the injection surface
			bundlespec.clear();
			BOOST_FOREACH(Track &track, tracks) {
				std::pair<double, double> steps =
				    surface->GetIntersection(track.GetPos(), track.GetDir());
				if (!std::isfinite(steps.first))
					continue;
				track.advance(steps.first);
				if (track.energy > 0)
					bundlespec.push_back(std::make_pair(track.radius(*primary), track.energy));
			}
			
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
