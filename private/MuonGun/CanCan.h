
#ifndef MUONGUN_CANCAN_H_INCLUDED
#define MUONGUN_CANCAN_H_INCLUDED

#include <MuonGun/Distribution.h>
#include <MuonGun/SingleMuonFlux.h>
#include <MuonGun/MultiplicityFraction.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <dataclasses/physics/I3Particle.h>
#include <MuonGun/I3MuonGun.h>

#include <boost/tuple/tuple.hpp>

namespace I3MuonGun {

class BundleInjector : public Distribution {
public:

	void SetRandomService(I3RandomServicePtr p);
	
	void SetSurface(CylinderPtr p);
	CylinderConstPtr GetSurface() const { return surface_; }
	
	void SetFlux(SingleMuonFluxPtr p);
	SingleMuonFluxConstPtr GetFlux() const { return flux_; }
	
	void SetMultiplicity(MultiplicityFractionPtr p);
	MultiplicityFractionConstPtr GetMultiplicity() const { return multiplicity_; }
	
	double GetTotalRate() const;
	
	// Draw a sample from the distribution of shower impact points,
	// returning the shower core, multiplicity, and rate contribution 
	boost::tuple<I3Particle, unsigned, double> GenerateAxis() const;
private:
	void CalculateMaxFlux();
	
	CylinderPtr surface_;
	SingleMuonFluxPtr flux_;
	MultiplicityFractionPtr multiplicity_;
	
	double maxFlux_;
	mutable double totalRate_;
};

}

#endif

