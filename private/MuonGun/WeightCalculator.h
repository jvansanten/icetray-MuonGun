
#ifndef I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED
#define I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <icetray/I3PointerTypedefs.h>

namespace I3MuonGun {

I3_FORWARD_DECLARATION(Generator);
I3_FORWARD_DECLARATION(SamplingSurface);
I3_FORWARD_DECLARATION(Flux);
I3_FORWARD_DECLARATION(RadialDistribution);
I3_FORWARD_DECLARATION(EnergyDistribution);

struct BundleModel {
	// BundleModel() {}
	BundleModel(FluxConstPtr f, RadialDistributionConstPtr r, EnergyDistributionConstPtr e)
	    : flux(f), radius(r), energy(e) {}
	FluxConstPtr flux;
	RadialDistributionConstPtr radius;
	EnergyDistributionConstPtr energy;
};

/**
 * @brief Utility class to calculate weights for muon bundles
 */
class WeightCalculator {
public:
	/**
	 * @param[in] s      Surface at which to calculate the weights
	 * @param[in] flux   Total flux in the desired model
	 * @param[in] radius Radial distribution in the desired model
	 * @param[in] energy Energy distribution in the desired model
	 * @param[in] g      Generation scheme according by which the
	 *                   events were generated.
	 */
	WeightCalculator(SamplingSurfacePtr s, FluxPtr flux,
	    RadialDistributionPtr radius, EnergyDistributionPtr energy, GenerationProbabilityPtr g)
	    : surface_(s), flux_(flux), radius_(radius), energy_(energy), generator_(g) {}
	
	/**
	 * Calculate a weight corresponding to the frequency with which
	 * the given bundle configuration appears in the flux model.
	 *
	 * @param[in] axis   The bundle axis
	 * @param[in] bundle The radial offset and energy of each muon
	 *                   in the bundle
	 * @returns a weight in units of @f$ [s^{-1}] @f$
	 */
	double GetWeight(const I3Particle &axis, const BundleConfiguration &bundle) const;
private:
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	RadialDistributionPtr radius_;
	EnergyDistributionPtr energy_;
	GenerationProbabilityPtr generator_;
};

}

#endif // I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED