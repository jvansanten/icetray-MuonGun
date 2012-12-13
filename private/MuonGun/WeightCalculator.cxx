
#include <MuonGun/WeightCalculator.h>
#include <MuonGun/I3MuonGun.h>
#include <boost/foreach.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
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
		
		// Harvest the direct daughters of the primary. This method is
		// appropriate only for the output of StaticCylinderInjector
		BundleConfiguration bundlespec;
		BOOST_FOREACH(const I3Particle &p, std::make_pair(mctree->begin(primary), mctree->end(primary))) {
			if (p.GetType() != I3Particle::MuMinus)
				continue;
			bundlespec.push_back(std::make_pair(p.GetPos().CalcDistance(primary->GetPos()), p.GetEnergy()));
		}
		
		double weight = calculator_->GetWeight(*primary, bundlespec);
		
		frame->Put("Weight", boost::make_shared<I3Double>(weight));
		PushFrame(frame);
	}
private:
	boost::shared_ptr<WeightCalculator> calculator_;
};

}

I3_MODULE(I3MuonGun::WeightCalculatorModule);
