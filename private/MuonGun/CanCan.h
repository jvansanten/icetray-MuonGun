
#ifndef MUONGUN_CANCAN_H_INCLUDED
#define MUONGUN_CANCAN_H_INCLUDED

#include <MuonGun/RadialDistribution.h>
#include <dataclasses/physics/I3Particle.h>
#include <MuonGun/I3MuonGun.h>

namespace I3MuonGun {

class SingleMuonFlux : public Distribution {
public:
	virtual ~SingleMuonFlux();	
	virtual double operator()(double depth, double cos_theta) const = 0;
};

I3_POINTER_TYPEDEFS(SingleMuonFlux);

// From Becherini et al.
class BMSSSingleMuonFlux : public SingleMuonFlux {
public:
	BMSSSingleMuonFlux();
	double operator()(double depth, double cos_theta) const;
private:
	double k0a_, k0b_, k1a_, k1b_;
};

// Becherini et al, re-fit to CORSIKA + MMC simulation
class AdHocSingleMuonFlux : public SingleMuonFlux {
public:
	AdHocSingleMuonFlux();
	double operator()(double depth, double cos_theta) const;
private:
	double a0_, a1_, b0_, b1_, c0_, c1_, d0_, d1_;
};

class MultiplicityFraction : public Distribution {
public:	
	double operator()(double depth, double cos_theta, unsigned multiplicity) const;
private:
	
};

I3_POINTER_TYPEDEFS(MultiplicityFraction);

class EnergyDistribution : public Distribution {
public:
	virtual ~EnergyDistribution() {};
	virtual double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius, double energy) const = 0;
	virtual Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(EnergyDistribution);

class PotemkinEnergyDistribution : public EnergyDistribution {
public:
	double operator()(double depth, double cos_theta, 
	    unsigned multiplicity, double radius, double energy) const;
	Sample Generate(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
};

class CanCan : public Distribution {
public:
	CanCan(double radius, double height, unsigned maxMultiplicity=1);
	
	void SetRandomService(I3RandomServicePtr r);
	
	// Fill sampled muons into the given vector, and return
	// the probability of obtaining the sampled configuration
	double GenerateBundle(std::vector<I3Particle> &) const;
	
	// Sample a shower axis and multiplicity
	std::pair<I3Particle, unsigned> GenerateAxis() const;
	
	// Integrate the flux parametrization over the sampling
	// surface to obtain the total rate in Hz
	double GetTotalRate() const;
	
	// Translate an IceCube z-coordinate into a depth in kmwe
	double GetDepth(double z) const;
	
private:
	I3Position GetImpactPoint(const I3Direction &dir) const;
	double GetDifferentialArea(double cos_zenith) const;
	double GetDifferentialTopArea(double cos_zenith) const;
	double GetDifferentialSideArea(double cos_zenith) const;
	
	double CalculateTotalRate(unsigned multiplicity) const;
	
	I3MuonGun::Cylinder surface_;
	SingleMuonFluxPtr flux_;
	MultiplicityFractionPtr multiFraction_;
	RadialDistributionPtr radialDistribution_;
	EnergyDistributionPtr energyDistribution_;
	double density_;
	double radius_, height_, originDepth_;
	unsigned maxMultiplicity_;
	double maxFlux_;
	mutable double totalRate_;
};

}

#endif

