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
 * sampled only on a small inner surface, while the surface scales up to full size
 * for potentially bright muons. This technique can be used to efficiently simulate
 * background for an event selection that requires a thick veto for dim events
 * (where the rates are also highest) but becomes more accepting for bright events.
 */
class EnergyDependentSurfaceInjector : public Generator {
public:
	EnergyDependentSurfaceInjector();
	
	// GenerationProbability interface
	SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const;
	double GetGenerationProbability(double depth, double coszen, const BundleConfiguration &bundle) const;
	GenerationProbabilityPtr Clone() const;
	
	// Generator interface
	void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	
	/** 
	 * Scale the sampling cylinder to a size appropriate for
	 * the given maximum muon energy
	 * by the given energy. This is not necessarily the fastest.
	 */
	SamplingSurfaceConstPtr GetSurface(double energy) const;
	/** 
	 * Integrate the flux to get the total rate on the surface determined
	 * by the given energy. This is not necessarily the fastest.
	 */
	double GetTotalRate(double energy) const;
private:
	CylinderPtr surface_;
	FluxPtr flux_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	RadialDistributionPtr radialDistribution_;
};

}

#endif
