/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_ENERGYDEPENDENTSURFACEINJECTOR_H_INCLUDED
#define I3MUONGUN_ENERGYDEPENDENTSURFACEINJECTOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <boost/function.hpp>

namespace I3MuonGun {

I3_FORWARD_DECLARATION(Cylinder);
I3_FORWARD_DECLARATION(SamplingSurface);
I3_FORWARD_DECLARATION(Flux);
I3_FORWARD_DECLARATION(OffsetPowerLaw);
I3_FORWARD_DECLARATION(RadialDistribution);

/**
 * @brief A rejection-sampling Generator with energy-dependent sampling surface
 *
 * EnergyDependentSurfaceInjector samples bundle impact points, angles, multiplicities,
 * and radial distributions at their natural frequencies, but scales the sampling
 * surface based on the highest-energy muon in the bundle: dim, low-energy muons are
 * targeted only at a small inner surface, while the surface scales up to full size
 * for potentially bright muons. This technique can be used to efficiently simulate
 * background for an event selection that requires a thick veto for dim events
 * (where the rates are also highest) but becomes more accepting for bright events.
 */
class EnergyDependentSurfaceInjector : public Generator {
public:
	EnergyDependentSurfaceInjector(FluxPtr flux=FluxPtr(), RadialDistributionPtr radius=RadialDistributionPtr(),
	    boost::shared_ptr<OffsetPowerLaw> energies=boost::shared_ptr<OffsetPowerLaw>(),
	    boost::function<SamplingSurfacePtr (double)> scaling=boost::function<SamplingSurfacePtr (double)>());
	
	// GenerationProbability interface
	SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const;
	double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	GenerationProbabilityPtr Clone() const;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const;
	
	// Generator interface
	void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	
	boost::function<SamplingSurfacePtr (double)> GetScaling() const { return scalingFunction_; }
	void SetScaling(boost::function<SamplingSurfacePtr (double)> f) { scalingFunction_ = f; }
	
	FluxConstPtr GetFlux() const { return flux_; }
	void SetFlux(FluxPtr f) { flux_ = f; }
	
	boost::shared_ptr<const OffsetPowerLaw> GetEnergyDistribution() const { return energyGenerator_; }
	void SetEnergyDistribution(boost::shared_ptr<OffsetPowerLaw> f) { energyGenerator_ = f; }
	
	RadialDistributionConstPtr GetRadialDistribution() const { return radialDistribution_; }
	void SetRadialDistribution(RadialDistributionPtr f) { radialDistribution_ = f; }
	
	/** 
	 * Scale the sampling cylinder to a size appropriate for
	 * the given maximum muon energy
	 * by the given energy. This is not necessarily the fastest.
	 */
	SamplingSurfacePtr GetTargetSurface(double energy) const;
	/** 
	 * Integrate the flux to get the total rate on the surface.
	 * This is not necessarily the fastest.
	 */
	double GetTotalRate(SamplingSurfaceConstPtr surface) const;
private:
	boost::function<SamplingSurfacePtr (double)> scalingFunction_;
	SamplingSurfacePtr injectionSurface_;
	FluxPtr flux_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	RadialDistributionPtr radialDistribution_;
};

}

#endif
