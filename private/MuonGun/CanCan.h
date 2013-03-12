/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef MUONGUN_CANCAN_H_INCLUDED
#define MUONGUN_CANCAN_H_INCLUDED

#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <dataclasses/physics/I3Particle.h>

#include <boost/tuple/tuple.hpp>

namespace I3MuonGun {

/**
 * @brief A simple rejection-sampling Generator
 *
 * StaticSurfaceInjector samples bundle impact points, angles, multiplicities,
 * and radial distributions at their natural frequencies on a fixed surface
 * using a brain-dead acceptance/rejection technique. Energies, on the other hand,
 * are sampled from an OffsetPowerLaw, which both boosts efficiency and allows
 * a measure of control over the generated energy distribution.
 */
class StaticSurfaceInjector : public Generator {
public:
	StaticSurfaceInjector();
	StaticSurfaceInjector(SamplingSurfacePtr surface, FluxPtr flux,
	    boost::shared_ptr<OffsetPowerLaw> edist, RadialDistributionPtr rdist);
	
	// Generator Interface
	void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	GenerationProbabilityPtr Clone() const;
	double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const { return surface_; }
	
	void SetSurface(SamplingSurfacePtr p);
	SamplingSurfacePtr GetSurface() { return surface_; }
	
	void SetFlux(FluxPtr p);
	FluxPtr GetFlux() { return flux_; }
	
	void SetRadialDistribution(RadialDistributionPtr r) { radialDistribution_ = r; }
	RadialDistributionPtr GetRadialDistribution() { return radialDistribution_; }
	
	void SetEnergyDistribution(boost::shared_ptr<OffsetPowerLaw> e) { energyGenerator_ = e; }
	boost::shared_ptr<OffsetPowerLaw> GetEnergyDistribution() { return energyGenerator_; }
	
	/**
	 * Integrate the configured flux over the sampling surface, summing over
	 * all allowed multiplicities.
	 *
	 * @returns a rate in units of @f$ [s^{-1}] @f$
	 */
	double GetTotalRate() const;
private:
	/**
	 * Draw a sample from the distribution of shower impact points
	 *
	 * @returns shower axis and multiplicity
	 */
	 void GenerateAxis(I3RandomService &rng, std::pair<I3Particle, unsigned> &axis) const;
	/**
	 * Distribute the given number of muons in the transverse plane
	 * and draw an energy for each
	 */
	void FillMCTree(I3RandomService &rng, const std::pair<I3Particle, unsigned> &axis, I3MCTree &, BundleConfiguration &) const;
	
	void CalculateMaxFlux();
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	RadialDistributionPtr radialDistribution_;
	
	double maxFlux_;
	mutable double totalRate_;

};

}

#endif

