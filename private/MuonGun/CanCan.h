
#ifndef MUONGUN_CANCAN_H_INCLUDED
#define MUONGUN_CANCAN_H_INCLUDED

#include <MuonGun/Distribution.h>
#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <dataclasses/physics/I3Particle.h>

#include <boost/tuple/tuple.hpp>

namespace I3MuonGun {

class StaticSurfaceInjector : public Generator {
public:
	StaticSurfaceInjector();
	
	// Generator Interface
	void Generate(I3MCTree &tree, BundleConfiguration &bundle) const;
	GenerationProbabilityPtr Clone() const;
	double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	SurfaceConstPtr GetInjectionSurface() const { return surface_; }
	
	void SetRandomService(I3RandomServicePtr p);
	
	void SetSurface(SamplingSurfacePtr p);
	SamplingSurfacePtr GetSurface() { return surface_; }
	
	void SetFlux(FluxPtr p);
	FluxPtr GetFlux() { return flux_; }
	
	void SetRadialDistribution(RadialDistributionPtr r) { radialDistribution_ = r; }
	RadialDistributionPtr GetRadialDistribution() { return radialDistribution_; }
	
	void SetEnergyDistribution(boost::shared_ptr<OffsetPowerLaw> e) { energyGenerator_ = e; }
	boost::shared_ptr<OffsetPowerLaw> GetEnergyDistribution() { return energyGenerator_; }
	
	double GetTotalRate() const;
private:
	// Draw a sample from the distribution of shower impact points,
	// returning the shower core and multiplicity
	void GenerateAxis(std::pair<I3Particle, unsigned> &axis) const;
	// Generate and distribute the given number of muons over the transverse plane
	void FillMCTree(const std::pair<I3Particle, unsigned> &axis, I3MCTree &, BundleConfiguration &) const;
	
	void CalculateMaxFlux();
	
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	RadialDistributionPtr radialDistribution_;
	
	double maxFlux_;
	mutable double totalRate_;

};

class BundleInjector : public Distribution {
public:

	void SetRandomService(I3RandomServicePtr p);
	
	void SetSurface(CylinderPtr p);
	CylinderConstPtr GetSurface() const { return surface_; }
	
	void SetFlux(FluxPtr p);
	FluxConstPtr GetFlux() const { return flux_; }
	
	double GetTotalRate() const;
	
	// Draw a sample from the distribution of shower impact points,
	// returning the shower core, multiplicity, and rate contribution 
	boost::tuple<I3Particle, unsigned, double> GenerateAxis() const;
private:
	void CalculateMaxFlux();
	
	CylinderPtr surface_;
	FluxPtr flux_;
	
	double maxFlux_;
	mutable double totalRate_;
};

}

#endif

