
#include "MuonGun/MuonPropagator.h"

#include "PROPOSAL/Propagate.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/MathModel.h"

namespace I3MuonGun {

MuonPropagator::MuonPropagator(const std::string &medium, double ecut, double vcut, double rho)
	: propagator_(new Propagate(medium, ecut, vcut, "mu", rho))
{
	propagator_->sdec      = true; // stopped muon decay
	propagator_->contiCorr = true; // continuous randomization
	propagator_->exactTime = true; // exact local time
	propagator_->molieScat = true; // Moliere scattering
	
	// LPM suppression
	propagator_->get_cros()->set_lpm(true); 
	// Kelner, Kokoulin, and Petrukhin parametrization
	propagator_->get_cros()->get_bremsstrahlung()->set_form(1);
	// Abramowicz Levin Levy Maor parametrization
	propagator_->get_cros()->get_photonuclear()->set_form(3);
	// ALLM 97 (rather than 91)
	propagator_->get_cros()->get_photonuclear()->set_bb(2);
	// Butkevich- Mikhailov nuclear structure function
	propagator_->get_cros()->get_photonuclear()->set_shadow(2);

	std::ostringstream prefix;
	prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";
	propagator_->interpolate("all", prefix.str());
}

MuonPropagator::~MuonPropagator()
{
	delete propagator_;
}

void
MuonPropagator::SetSeed(int seed)
{
	MathModel::set_seed(seed);
}

std::string
MuonPropagator::GetName(const I3Particle &p)
{
        std::string name;
	
	switch (p.GetType()) {
		case I3Particle::MuMinus:
			name="mu-";
			break;
		case I3Particle::MuPlus:
			name="mu+";
			break;
		default:
			break;
	}
	
	return name;
}

I3Particle
MuonPropagator::propagate(const I3Particle &p, double distance)
{
	I3Particle endpoint(p);
	
	// propagator_.get_output()->DEBUG=true;
	propagator_->get_output()->initDefault(0, 0, GetName(p), p.GetTime()/I3Units::second,
	    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
	    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);
	
	PROPOSALParticle *pp = propagator_->get_particle();
	if (propagator_->propagateTo(distance/I3Units::cm, p.GetEnergy()/I3Units::MeV) > 0)
		endpoint.SetEnergy(pp->e*I3Units::MeV);
	else
		endpoint.SetEnergy(0);
	
	endpoint.SetPos(pp->x*I3Units::cm, pp->y*I3Units::cm, pp->z*I3Units::cm);
	endpoint.SetThetaPhi(pp->theta*I3Units::degree, pp->phi*I3Units::degree);
	endpoint.SetLength(pp->r*I3Units::cm);
	endpoint.SetTime(pp->t*I3Units::second);
	
	return endpoint;
}

}