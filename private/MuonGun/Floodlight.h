
#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/EnergyDistribution.h>

namespace I3MuonGun {

class Floodlight : public Generator {
public:
	Floodlight();
	
	// Generator Interface
	void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	GenerationProbabilityPtr Clone() const;
	double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const { return surface_; }
	
private:
	SamplingSurfacePtr surface_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;

};

}