#ifndef MUONGUN_RADIALDISTRIBUTION_H
#define MUONGUN_RADIALDISTRIBUTION_H

#include <icetray/I3FrameObject.h>
#include <phys-services/I3RandomService.h>

class I3Position;

namespace I3MuonGun {
	
// These classes are effectively a re-implemenation of MUPAGE:
// http://arxiv.org/abs/0802.0562v2
	
// Get the depth of the given IceCube coordinate in kilometers water-equivalent
// double GetDepth(const I3Position &);

class Distribution : public I3FrameObject {
public:
	virtual ~Distribution() {};
	I3RandomServicePtr GetRandomService() const { return rng_; }
	void SetRandomService(I3RandomServicePtr r) { rng_ = r; }
	
protected:
	I3RandomServicePtr rng_;
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Distribution);

class RadialDistribution : public Distribution {
public:
	virtual ~RadialDistribution() {};
	virtual double GetGenerationProbability(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
	virtual double Generate(double depth, double cos_theta,
	    unsigned multiplicity) const = 0;
private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(RadialDistribution);

class BMSSRadialDistribution : public RadialDistribution {
public:
	BMSSRadialDistribution();
	virtual ~BMSSRadialDistribution() {};
	virtual double GetGenerationProbability(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	virtual double Generate(double depth, double cos_theta,
	    unsigned multiplicity) const;
private:
	double GetMeanRadius(double, double, unsigned) const;
	double GetShapeParameter(double, double, unsigned) const;
	double GetGenerationProbability(double, double, double) const;
	double rho0a_, rho0b_, rho1_, theta0_, f_, alpha0a_, alpha0b_, alpha1a_, alpha1b_;
	double rmax_;
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(BMSSRadialDistribution);

}

#endif // MUONGUN_RADIALDISTRIBUTION_H
