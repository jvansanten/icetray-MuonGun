
#include <MuonGun/Generator.h>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <phys-services/I3RandomService.h>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

namespace I3MuonGun {

GenerationProbability::~GenerationProbability() {}
Generator::~Generator() {}

double
GenerationProbability::GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return numEvents_*GetGenerationProbability(axis, bundle);
}

GenerationProbabilityCollection::GenerationProbabilityCollection(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{
	push_back(p1);
	push_back(p2);
}

double
GenerationProbabilityCollection::GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	double prob = 0.;
	BOOST_FOREACH(const value_type &p, *this)
		if (p)
			prob += p->GetGeneratedEvents(axis, bundle);
	
	return prob;
}

GenerationProbabilityPtr
GenerationProbabilityCollection::Clone() const
{
	return boost::make_shared<GenerationProbabilityCollection>(*this);
}

GenerationProbabilityPtr
operator*=(GenerationProbabilityPtr p, size_t n)
{
	p->SetTotalEvents(p->GetTotalEvents()*n);
	
	return p;
}

GenerationProbabilityPtr
operator*(GenerationProbabilityPtr p, size_t n)
{
	GenerationProbabilityPtr pn = p->Clone();
	pn *= n;
	
	return pn;
}

GenerationProbabilityPtr
operator*(size_t n, GenerationProbabilityPtr p)
{
	return p*n;
}

GenerationProbabilityPtr
operator+(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{
	boost::shared_ptr<GenerationProbabilityCollection> c1 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p1);
	boost::shared_ptr<GenerationProbabilityCollection> c2 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p2);
	if (c1 && c2) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		std::copy(c2->begin(), c2->end(), std::back_inserter(*c1));
		return c1;
	} else if (c1) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		c1->push_back(p2);
		return c1;
	} else if (c2) {
		c2 = boost::make_shared<GenerationProbabilityCollection>(*c2);
		c2->push_back(c1);
		return c2;
	}
	
	return boost::make_shared<GenerationProbabilityCollection>(p1, p2);
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
		
		rng_ = context_.Get<I3RandomServicePtr>();
		if (!rng_)
			log_fatal("No RandomService configured!");
		maxEvents_ = generator_->GetTotalEvents();
	}
	
	void DAQ(I3FramePtr frame)
	{
		I3MCTreePtr mctree = boost::make_shared<I3MCTree>();
		BundleConfiguration bundlespec;
		
		generator_->Generate(*rng_, *mctree, bundlespec);
		frame->Put(mctreeName_, mctree);
		
		PushFrame(frame);
		if (++numEvents_ >= maxEvents_)
			RequestSuspension();
	}
private:
	GeneratorPtr generator_;
	I3RandomServicePtr rng_;
	size_t maxEvents_, numEvents_;
	std::string mctreeName_;
};

}

I3_MODULE(I3MuonGun::GeneratorModule);
