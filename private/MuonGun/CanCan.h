
#include <MuonGun/RadialDistribution.h>
#include <dataclasses/physics/I3Particle.h>

namespace I3MuonGun {

class SingleMuonFlux : public Distribution {
public:
	SingleMuonFlux();
	
	double GetFlux(double depth, double cos_theta) const;
private:
	double a0_, a1_, b0_, b1_, c0_, c1_, d0_, d1_;
};

class CanCan : public Distribution {
public:
	CanCan(double radius, double height);
	
	I3Particle Generate() const;
	double GetLivetime() const;
private:
	I3Position GetImpactPoint(const I3Direction &dir) const;
	double GetDifferentialArea(double cos_zenith) const;
	double GetDifferentialTopArea(double cos_zenith) const;
	double GetDifferentialSideArea(double cos_zenith) const;
	
	double GetDepth(double z) const;
	
	SingleMuonFlux flux_;
	double density_;
	double radius_, height_, originDepth_;
	double maxFlux_;
};

}
