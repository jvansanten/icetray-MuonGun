
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