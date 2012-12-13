
#include <MuonGun/Generator.h>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <phys-services/I3RandomService.h>
#include <boost/make_shared.hpp>

namespace I3MuonGun {

Generator::~Generator() {}

double
Generator::GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return numEvents_*GetGenerationProbability(axis, bundle);
}

class GeneratorModule : public I3Module {
public:
	GeneratorModule(const I3Context &ctx) : I3Module(ctx), maxEvents_(0), numEvents_(0)
	{
		AddOutBox("OutBox");
		AddParameter("Generator", "", generator_);
		
		mctreeName_ = "I3MCTree";
	}
	
	void Configure()
	{
		GetParameter("Generator", generator_);
		
		generator_->SetRandomService(context_.Get<I3RandomServicePtr>());
		maxEvents_ = generator_->GetTotalEvents();
	}
	
	void DAQ(I3FramePtr frame)
	{
		I3MCTreePtr mctree = boost::make_shared<I3MCTree>();
		BundleConfiguration bundlespec;
		
		generator_->Generate(*mctree, bundlespec);
		// TODO: propagate muons
		frame->Put(mctreeName_, mctree);
		
		PushFrame(frame);
		if (++numEvents_ >= maxEvents_)
			RequestSuspension();
	}
private:
	GeneratorPtr generator_;
	size_t maxEvents_, numEvents_;
	std::string mctreeName_;
};

}

I3_MODULE(I3MuonGun::GeneratorModule);
