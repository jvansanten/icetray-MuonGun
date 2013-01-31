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
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return 0.;
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundlespec.size();
	
	double norm = generator_->GetGeneratedEvents(axis, bundlespec);
	double rate = (*flux_)(h, coszen, m)*surface_->GetDifferentialArea(coszen)/norm;
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

MuonBundleConverter::MuonBundleConverter(size_t maxMultiplicity, SamplingSurfacePtr surface)
    : maxMultiplicity_(maxMultiplicity),
    surface_(surface ? surface : boost::make_shared<Cylinder>(1600, 800))
{}

I3TableRowDescriptionPtr
MuonBundleConverter::CreateDescription(const I3MCTree&)
{
	I3TableRowDescriptionPtr desc(new I3TableRowDescription());
	
	desc->AddField<uint32_t>("multiplicity", "", "Number of muons in the bundle");
	desc->AddField<float>("depth", "km", "Vertical depth of intersection with the sampling surface");
	desc->AddField<float>("cos_theta", "", "Cosine of the shower zenith angle");
	desc->AddField<float>("energy", "GeV", "Muon energy at sampling surface",
	    maxMultiplicity_);
	desc->AddField<float>("radius", "m", "Perpendicular distance from of track "
	    "from the bundle axis at the sampling surface", maxMultiplicity_);
	
	return desc;
}

size_t
MuonBundleConverter::FillRows(const I3MCTree &mctree, I3TableRowPtr rows)
{
	I3MMCTrackListConstPtr mmctracks = currentFrame_->Get<I3MMCTrackListConstPtr>("MMCTrackList");
	if (!mmctracks)
		log_fatal("I3MMCTrackList missing!");
	
	const I3MCTree::iterator primary = mctree.begin();
	std::pair<double, double> steps =
	    surface_->GetIntersection(primary->GetPos(), primary->GetDir());
	if (steps.first > 0) {
		rows->Set<float>("depth", GetDepth(primary->GetPos().GetZ() + steps.first*primary->GetDir().GetZ()));
		rows->Set<float>("cos_theta", cos(primary->GetDir().GetZenith()));
	}
	
	uint32_t m = 0;
	float *energies = rows->GetPointer<float>("energy");
	float *radii = rows->GetPointer<float>("radius");
	
	BOOST_FOREACH(const Track &track, Track::Harvest(mctree, *mmctracks)) {
		std::pair<double, double> steps =
		    surface_->GetIntersection(track.GetPos(), track.GetDir());
		float energy = track.GetEnergy(steps.first);
		if (energy > 0) {
			if (m < maxMultiplicity_) {
				energies[m] = energy;
				radii[m] = GetRadius(*primary, track.GetPos(steps.first));
			}
			m++;
		}
	}
	
	rows->Set("multiplicity", m);
	
	return 1;
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
		
		// This should be a fixed surface
		SamplingSurfaceConstPtr surface = generator_->GetInjectionSurface(*primary, bundlespec);
		
		double rate = 0.;
		std::pair<double, double> steps =
		    surface->GetIntersection(primary->GetPos(), primary->GetDir());
		if (bundlespec.size() > 0 && std::isfinite(steps.first)) {
			double h = GetDepth(primary->GetPos().GetZ() + steps.first*primary->GetDir().GetZ());
			double coszen = cos(primary->GetDir().GetZenith());
			unsigned m = bundlespec.size();
			
			double norm = generator_->GetGeneratedEvents(*primary, bundlespec);
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
