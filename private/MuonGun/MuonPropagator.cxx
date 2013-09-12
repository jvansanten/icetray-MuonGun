/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

#include "MuonGun/MuonPropagator.h"

#include "PROPOSAL/Propagate.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/MathModel.h"

#include "PROPOSAL/Medium.h"
#include "PROPOSAL/BremsStochastic.h"
#include "PROPOSAL/EpairStochastic.h"
#include "PROPOSAL/PhotoStochastic.h"
#include "PROPOSAL/IonizStochastic.h"

namespace I3MuonGun {

MuonPropagator::MuonPropagator(const std::string &medium, double ecut, double vcut, double rho)
	: propagator_(new Propagate(medium, ecut, vcut, "mu", rho))
{
	
	propagator_->sdec      = true; // stopped muon decay
	propagator_->exactTime = true; // exact local time
	propagator_->molieScat = true; // Moliere scattering
	
	// Turn on continuous randomization if no absolute
	// energy cutoff set
	propagator_->contiCorr = ecut < 0;
	propagator_->contiCorr = false;
	
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

inline std::string
GetMMCName(I3Particle::ParticleType pt)
{
	std::string name;
	
	switch (pt) {
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

std::string
MuonPropagator::GetName(const I3Particle &p)
{
	return GetMMCName(p.GetType());
}

typedef boost::bimaps::bimap<boost::bimaps::multiset_of<I3Particle::ParticleType>, boost::bimaps::multiset_of<int> > particle_type_conversion_t;

static const particle_type_conversion_t fromRDMCTable =
boost::assign::list_of<particle_type_conversion_t::relation>
(I3Particle::unknown, -100)
(I3Particle::Gamma, 1)
(I3Particle::EPlus, 2)
(I3Particle::EMinus, 3)
(I3Particle::Nu, 4)
(I3Particle::MuPlus, 5)
(I3Particle::MuMinus, 6)
(I3Particle::Pi0, 7)
(I3Particle::PiPlus, 8)
(I3Particle::PiMinus, 9)
(I3Particle::KPlus, 11)
(I3Particle::KMinus, 12)
(I3Particle::PPlus, 14)
(I3Particle::PMinus, 15)
(I3Particle::TauPlus, 33)
(I3Particle::TauMinus, 34)
(I3Particle::Monopole, 41)
(I3Particle::NuE, 201)
(I3Particle::NuMu, 202)
(I3Particle::NuTau, 203)
(I3Particle::NuEBar, 204)
(I3Particle::NuMuBar, 205)
(I3Particle::NuTauBar, 206)
(I3Particle::Brems, 1001)
(I3Particle::DeltaE, 1002)
(I3Particle::PairProd, 1003)
(I3Particle::NuclInt, 1004)
(I3Particle::MuPair, 1005)
(I3Particle::Hadrons, 1006);

inline I3Particle
to_I3Particle(const PROPOSALParticle *pp)
{
	I3Particle p;
	particle_type_conversion_t::right_const_iterator it =
	    fromRDMCTable.right.find(abs(pp->type));
	if (it == fromRDMCTable.right.end())
		log_fatal("unknown RDMC code \"%i\" cannot be converted to a I3Particle::ParticleType.", pp->type);
	else
		p.SetType(it->second);
	p.SetLocationType(I3Particle::InIce);
	p.SetPos(pp->x*I3Units::cm, pp->y*I3Units::cm, pp->z*I3Units::cm);
	p.SetTime(pp->t*I3Units::s);
	p.SetThetaPhi(pp->theta*I3Units::deg, pp->phi*I3Units::deg);
	p.SetLength(pp->l*I3Units::cm);
	p.SetEnergy(pp->e*I3Units::MeV);
	
	return p;
}

/** Differential stochastic rate: d^2N/dv/dx [1/m] */
double
MuonPropagator::GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type) const
{
	propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
	// Check kinematics
	if (fraction <= 0 || energy*(1-fraction) <= propagator_->get_particle()->m*I3Units::MeV)
		return 0.;
	propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	propagator_->get_cros()->get_ionization()->setEnergy();
	
	double contrib;
	double rate = 0.;
	// Separate contributions from each element for brems/epair/photonuclear interactions
	for (int i=0; i < propagator_->get_cros()->get_medium()->get_numCompontents(); i++) {
		propagator_->get_cros()->set_component(i);
		if (std::isfinite(contrib = propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
		if (std::isfinite(contrib = propagator_->get_cros()->get_epairproduction()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
		if (std::isfinite(contrib = propagator_->get_cros()->get_photonuclear()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
	}
	// Only one bulk ionization contribution
	if (std::isfinite(contrib = propagator_->get_cros()->get_ionization()->get_Stochastic()->function(fraction)) && contrib > 0)
		rate += contrib;
	// printf("brems dN/dx: %e\n", propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx());
	// printf("epair dN/dx: %e\n", propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx());
	// printf("photo dN/dx: %e\n", propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx());
	// printf("ioniz dN/dx: %e\n", propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx());
	
	return rate*(I3Units::m/I3Units::cm);
}

/** total stochastic rate: dN/dx [1/m] */
double
MuonPropagator::GetTotalStochasticRate(double energy, I3Particle::ParticleType type) const
{
	propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
	propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	propagator_->get_cros()->get_ionization()->setEnergy();
	
	double rate = 0;
	rate += propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx();
	
	return rate*(I3Units::m/I3Units::cm);
}

I3Particle
MuonPropagator::propagate(const I3Particle &p, double distance, boost::shared_ptr<std::vector<I3Particle> > losses)
{
	I3Particle endpoint(p);
	
	// propagator_.get_output()->DEBUG=true;
	if (losses) {
		propagator_->get_output()->I3flag = true;
		propagator_->get_output()->initF2000(0, 0, GetName(p), p.GetTime()/I3Units::second,
		    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
		    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);
	} else
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
	
	if (losses) {
		std::vector<PROPOSALParticle*> &history = propagator_->get_output()->I3hist;
		BOOST_FOREACH(PROPOSALParticle *pp, history) {
			losses->push_back(to_I3Particle(pp));
			delete pp;
		}
		history.clear();
		propagator_->get_output()->I3flag = false;
	}
	
	return endpoint;
}

void Crust::AddLayer(boost::shared_ptr<Surface> s, boost::shared_ptr<MuonPropagator> p)
{
	boundaries_.push_back(s);
	propagators_.push_back(p);
}

I3Particle
Crust::Ingest(const I3Particle &p)
{
	I3Particle propped(p);
	double l = 0;
	for (unsigned i = 0; (propped.GetEnergy() > 0) && (i < boundaries_.size()); i++) {
		double dx = boundaries_[i]->GetIntersection(propped.GetPos(), propped.GetDir()).first;
		if (dx > 0)
			propped = (i > 0 ? propagators_[i-1] : defaultPropagator_)->propagate(propped, dx);
		// Force lengths to measure the distance back to the outermost surface
		if (i > 0)
			l += std::min(dx, propped.GetLength());
	}
	propped.SetLength(l);
	
	return propped;
}


}
