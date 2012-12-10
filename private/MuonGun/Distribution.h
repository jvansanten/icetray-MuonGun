
#ifndef MUONGUN_DISTRIBUTION_H
#define MUONGUN_DISTRIBUTION_H

#include <icetray/I3FrameObject.h>
#include <phys-services/I3RandomService.h>

namespace I3MuonGun {

class Distribution : public I3FrameObject {
public:
	virtual ~Distribution();
	I3RandomServicePtr GetRandomService() const { return rng_; }
	virtual void SetRandomService(I3RandomServicePtr r) { rng_ = r; }
	
protected:
	I3RandomServicePtr rng_;
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Distribution);

struct Sample {
	double value, prob;
	Sample(double v, double p) : value(v), prob(p) {}
};

}

#endif
