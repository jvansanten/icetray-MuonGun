
#ifndef MUONGUN_MUONPROPAGATOR_H_INCLUDED
#define MUONGUN_MUONPROPAGATOR_H_INCLUDED

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Particle.h"
#include "MuonGun/I3MuonGun.h"

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

// A set of nested media layers
class Crust {
public:
	Crust(boost::shared_ptr<MuonPropagator> defaultPropagator) : defaultPropagator_(defaultPropagator) {};
	
	// Add an inner layer
	void AddLayer(boost::shared_ptr<Surface>, boost::shared_ptr<MuonPropagator>);
	// Propagate a muon to the outer boundary of the innermost layer
	I3Particle Ingest(const I3Particle &p);
private:
	boost::shared_ptr<MuonPropagator> defaultPropagator_;
	std::vector<boost::shared_ptr<Surface> > boundaries_;
	std::vector<boost::shared_ptr<MuonPropagator> > propagators_;
};

}

#endif // MUONGUN_MUONPROPAGATOR_H_INCLUDED