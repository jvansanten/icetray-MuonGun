/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 *$Revision$
 * $Date$
 */

#ifndef I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED
#define I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED

#include <MuonGun/Generator.h>
#include <icetray/I3PointerTypedefs.h>

// only for GetMuonsAtSurface
#include <icetray/I3Frame.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>

#include <tableio/I3Converter.h>

namespace I3MuonGun {

I3_FORWARD_DECLARATION(Generator);
I3_FORWARD_DECLARATION(SamplingSurface);
I3_FORWARD_DECLARATION(Flux);
I3_FORWARD_DECLARATION(RadialDistribution);
I3_FORWARD_DECLARATION(EnergyDistribution);

std::vector<I3Particle>
GetMuonsAtSurface(I3FramePtr frame, SurfaceConstPtr surface);

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
	 
	SamplingSurfacePtr GetSurface() { return surface_; }
	void SetSurface(SamplingSurfacePtr s) { surface_ = s; }
private:
	SamplingSurfacePtr surface_;
	FluxPtr flux_;
	RadialDistributionPtr radius_;
	EnergyDistributionPtr energy_;
	GenerationProbabilityPtr generator_;
};

class MuonBundleConverter : public I3ConverterImplementation<I3MCTree> {
public:
	MuonBundleConverter(size_t maxMultiplicity=25, SamplingSurfacePtr surface=SamplingSurfacePtr());
	I3TableRowDescriptionPtr CreateDescription(const I3MCTree&);
	size_t FillRows(const I3MCTree&, I3TableRowPtr);
private:
	size_t maxMultiplicity_;
	SamplingSurfacePtr surface_;
};

}

#endif // I3MUONGUN_WEIGHTCALCULATOR_H_INCLUDED