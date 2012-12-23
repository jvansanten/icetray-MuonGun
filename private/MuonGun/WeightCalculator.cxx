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

namespace I3MuonGun {

double
WeightCalculator::GetWeight(const I3Particle &axis, const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect sampling surface. Bail.
	// if (steps.first == 0)
	// 	return 0.;
	
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

/**
 * @brief Interface between WeightCalculator and IceTray
 *
 * WeightCalculatorModule handles the details of extracting energies and
 * radial offsets of muons from an I3MCTree. 
 */
class WeightCalculatorModule : public I3Module {
public:
	WeightCalculatorModule(const I3Context &ctx) : I3Module(ctx)
	{
		AddOutBox("OutBox");
		AddParameter("WeightCalculator", "", calculator_);
	}
	
	void Configure()
	{
		GetParameter("WeightCalculator", calculator_);
		if (!calculator_)
			log_fatal("No weight calculator configured!");
	}
	
	void DAQ(I3FramePtr frame)
	{
		I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
		I3MCTree::iterator primary = mctree->begin();
		

		BundleConfiguration bundlespec;
		if (primary->GetType() == I3Particle::PPlus || primary->IsNucleus()) {
			// For complete air-shower simulation, harvest radial offsets and
			// energies at the sampling surface from the MMCTrackList
			I3MMCTrackListConstPtr mmctracks = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");
			if (!mmctracks)
				log_fatal("This appears to be CORSIKA simulation, but I have no MMCTrackList!");
			BOOST_FOREACH(const I3MMCTrack &track, *mmctracks) {
				// Get the position where the track crosses the sampling surface,
				// in a coordinate system where the shower plane is the
				// X-Y plane and the shower axis is the z axis.
				I3Position impact(track.GetXi()-primary->GetPos().GetX(),
				                  track.GetYi()-primary->GetPos().GetY(),
				                  track.GetZi()-primary->GetPos().GetZ());
				impact.RotateZ(-primary->GetDir().GetAzimuth());
				impact.RotateY(-primary->GetDir().GetZenith());
				impact.RotateZ(primary->GetDir().GetAzimuth());
				
				bundlespec.push_back(std::make_pair(impact.GetRho(), track.GetEi()));
			}
		} else {
			// Harvest the direct daughters of the primary. This method is
			// appropriate only for the output of StaticCylinderInjector
			BOOST_FOREACH(const I3Particle &p,
			    std::make_pair(mctree->begin(primary), mctree->end(primary))) {
				if (p.GetType() != I3Particle::MuMinus)
					continue;
				bundlespec.push_back(std::make_pair(
				    p.GetPos().CalcDistance(primary->GetPos()), p.GetEnergy()));
			}
		}
		
		double weight = calculator_->GetWeight(*primary, bundlespec);
		
		frame->Put(GetName(), boost::make_shared<I3Double>(weight));
		PushFrame(frame);
	}
private:
	boost::shared_ptr<WeightCalculator> calculator_;
};

}

I3_MODULE(I3MuonGun::WeightCalculatorModule);
