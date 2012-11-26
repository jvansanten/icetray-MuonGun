
#ifndef MUONGUN_MUONPROPAGATOR_H_INCLUDED
#define MUONGUN_MUONPROPAGATOR_H_INCLUDED

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Particle.h"

class Propagate;

namespace I3MuonGun {

class MuonPropagator {
public:
	MuonPropagator(const std::string &medium, double ecut=-1, double vcut=-1, double rho=1.0);
	~MuonPropagator();
	I3Particle propagate(const I3Particle &p, double distance);
	
	static void SetSeed(int seed);
	static std::string GetName(const I3Particle &p);
	
private:
	Propagate *propagator_;
};

}

#endif // MUONGUN_MUONPROPAGATOR_H_INCLUDED