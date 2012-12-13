
#ifndef I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED
#define I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

namespace I3MuonGun {

struct BundleModel {
	// BundleModel() {}
	BundleModel(FluxConstPtr f, RadialDistributionConstPtr r, EnergyDistributionConstPtr e)
	    : flux(f), radius(r), energy(e) {}
	FluxConstPtr flux;
	RadialDistributionConstPtr radius;
	EnergyDistributionConstPtr energy;
};

class WeightCalculator {
public:
	WeightCalculator(SamplingSurfacePtr s, FluxPtr flux,
	    RadialDistributionPtr radius, EnergyDistributionPtr energy, GeneratorPtr g)
	    : surface_(s), flux_(flux), radius_(radius), energy_(energy), generator_(g) {}
	
	double GetWeight(const I3Particle &, const BundleConfiguration &) const;
private:
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	RadialDistributionPtr radius_;
	EnergyDistributionPtr energy_;
	GeneratorPtr generator_;
};

}

#endif // I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED