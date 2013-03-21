
#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>
#include <MuonGun/EnergyDistribution.h>

namespace I3MuonGun {

class Floodlight : public Generator {
public:
	Floodlight();
	
	// Generator Interface
	virtual void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const;
	virtual GenerationProbabilityPtr Clone() const;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const;
	virtual double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
	virtual SamplingSurfaceConstPtr GetInjectionSurface() const { return surface_; }
	
private:
	double GetTotalRate() const;
	SamplingSurfacePtr surface_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	mutable double totalRate_;

};

}
