#ifndef I3MUONGUN_RADIALDISTRIBUTION_H
#define I3MUONGUN_RADIALDISTRIBUTION_H

#include <MuonGun/Distribution.h>
#include <photospline/I3SplineTable.h>

class I3Position;

namespace I3MuonGun {

class RadialDistribution : public Distribution {
public:
	virtual ~RadialDistribution() {};
	// TODO: short-circuit for multiplicity == 1
	virtual double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
	virtual Sample Generate(double depth, double cos_theta,
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
	virtual double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	virtual Sample Generate(double depth, double cos_theta,
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

class SplineRadialDistribution : public RadialDistribution, private I3SplineTable {
public:
	SplineRadialDistribution(const std::string&);
	double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity) const;
private:
	double RawProbability(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
};

}

#endif // I3MUONGUN_RADIALDISTRIBUTION_H
