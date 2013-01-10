/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/EnergyDistribution.h>
#include <phys-services/I3RandomService.h>

namespace I3MuonGun {

EnergyDistribution::~EnergyDistribution() {};

SplineEnergyDistribution::SplineEnergyDistribution(const std::string &singles, const std::string &bundles)
    : singles_(singles), bundles_(bundles)
{
	if (singles_.GetNDim() != 3u)
		log_fatal("'%s' does not appear to be a single-muon energy distribution", singles.c_str());
	if (bundles_.GetNDim() != 5u)
		log_fatal("'%s' does not appear to be a muon bundle energy distribution", bundles.c_str());
	// Extrapolate with a constant below the minimum supported energy
	minLogEnergy_ = std::min(singles_.GetExtents(2).first, bundles_.GetExtents(2).first);
}

double
SplineEnergyDistribution::operator()(double depth, double cos_theta, 
    unsigned multiplicity, double radius, double energy) const
{
	double coords[5] = {cos_theta, depth, multiplicity, radius, std::max(minLogEnergy_, std::log(energy))};
	double logprob;
	
	if (multiplicity < 2) {
		coords[2] = coords[4];
		if (singles_.Eval(coords, &logprob) != 0)
			return 0.;
	} else if (bundles_.Eval(coords, &logprob) != 0)
		return 0.;
	
	return std::exp(logprob);
}

double
SplineEnergyDistribution::Generate(I3RandomService &rng, double depth, double cos_theta, 
    unsigned multiplicity, double radius) const
{
	log_fatal("Sampling is not yet implemented");
	return 1.;
}

BMSSEnergyDistribution::BMSSEnergyDistribution() : 
    beta_(0.42), g0_(-0.232), g1_(3.961), e0a_(0.0304), e0b_(0.359), e1a_(-0.0077), e1b_(0.659),
    a0_(0.0033), a1_(0.0079), b0a_(0.0407), b0b_(0.0283), b1a_(-0.312), b1b_(6.124),
    q0_(0.0543), q1_(-0.365), c0a_(-0.069), c0b_(0.488), c1_(-0.117),
    d0a_(-0.398), d0b_(3.955), d1a_(0.012), d1b_(-0.35)
{}

double
BMSSEnergyDistribution::operator()(double depth, double cos_theta, 
    unsigned m, double r, double energy) const
{
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double bX = beta_*h/cos_theta;
	double g, eps;
	if (m == 1) {
		g = g0_*log(h) + g1_;
		eps = e0a_*exp(e0b_*h)/cos_theta + e1a_*h + e1b_;
	} else {
		m = std::min(m, 4u);
		double a = a0_*h + a1_;
		double b = (b0a_*m + b0b_)*h + (b1a_*m + b1b_);
		double q = q0_*h + q1_;
		g = a*r + b*(1 - 0.5*exp(q*r));
		double c = (c0a_*h + c0b_)*exp(c1_*r);
		double d = (d0a_*h + d0b_)*pow(r, d1a_*h + d1b_);
		eps = c*acos(cos_theta) + d;
	}
	double norm = (g-1)*pow(eps, g-1)*exp((g-1)*bX)*pow(1-exp(-bX), g-1);
	return norm*exp((1-g)*bX)*pow(energy + eps*(1-exp(-bX)), -g);
}

double
BMSSEnergyDistribution::Generate(I3RandomService &rng, double depth, double cos_theta, 
    unsigned multiplicity, double radius) const
{
	log_fatal("Sampling is not yet implemented");
	return 1.;
}

OffsetPowerLaw::OffsetPowerLaw(double gamma, double offset, double emin, double emax)
    : gamma_(gamma), offset_(offset), emin_(emin), emax_(emax)
{
	if (gamma <= 1)
		log_fatal("gamma must be > 1");
	nmin_ = std::pow(emin + offset, 1-gamma);
	nmax_ = std::pow(emax + offset, 1-gamma);
	norm_ = (1-gamma)/(nmax_ - nmin_);
}

double
OffsetPowerLaw::operator()(double energy) const
{
	if (energy <= emax_ && energy >= emin_)
		return norm_*std::pow(energy + offset_, -gamma_);
	else
		return 0.;
}

double
OffsetPowerLaw::Generate(I3RandomService &rng) const
{
	return std::pow(rng.Uniform()*(nmax_ - nmin_) + nmin_, 1./(1.-gamma_)) - offset_;
}

}