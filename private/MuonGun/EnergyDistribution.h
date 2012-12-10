
#ifndef I3MUONGUN_ENERGYDISTRIBUTION_H_INCLUDED
#define I3MUONGUN_ENERGYDISTRIBUTION_H_INCLUDED

#include <MuonGun/Distribution.h>
#include <icetray/I3Units.h>
#include <photospline/I3SplineTable.h>

namespace I3MuonGun {

class EnergyDistribution : public Distribution {
public:
	EnergyDistribution() : min_(I3Units::GeV), max_(I3Units::PeV) {}
	virtual ~EnergyDistribution();
	typedef double result_type;
	virtual double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius, double energy) const = 0;
	virtual Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
	
	double GetMax() const { return max_; }
	double GetMin() const { return min_; }
	void SetMax(double v) { max_ = v; }
	void SetMin(double v) { min_ = v; }
private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double min_, max_;
};

I3_POINTER_TYPEDEFS(EnergyDistribution);

class PotemkinEnergyDistribution : public EnergyDistribution {
public:
	double operator()(double depth, double cos_theta, 
	    unsigned multiplicity, double radius, double energy) const;
	Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
};

class SplineEnergyDistribution : public EnergyDistribution {
public:
	SplineEnergyDistribution(const std::string &singles, const std::string &bundles);
	double operator()(double depth, double cos_theta, 
	    unsigned multiplicity, double radius, double energy) const;
	Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
private:
	I3SplineTable singles_;
	I3SplineTable bundles_;
};

class OffsetPowerLaw : public Distribution {
public:
	OffsetPowerLaw(double gamma, double offset, double emin, double emax);
	typedef double result_type;
	double operator()(double energy) const;
	double Generate() const;
private:
	double gamma_, offset_;
	double emin_, emax_;
	double nmin_, nmax_, norm_;
};

}

#endif