
#include <MuonGun/EnergyDistribution.h>

namespace I3MuonGun {

EnergyDistribution::~EnergyDistribution() {};

double
PotemkinEnergyDistribution::operator()(double depth, double cos_theta, 
    unsigned multiplicity, double radius, double energy) const
{
	return 1.;
}

Sample
PotemkinEnergyDistribution::Generate(double depth, double cos_theta, 
    unsigned multiplicity, double radius) const
{
	return Sample(1., 1.);
}

SplineEnergyDistribution::SplineEnergyDistribution(const std::string &singles, const std::string &bundles)
    : singles_(singles), bundles_(bundles)
{
	if (singles_.GetNDim() != 3u)
		log_fatal("'%s' does not appear to be a single-muon energy distribution", singles.c_str());
	if (bundles_.GetNDim() != 5u)
		log_fatal("'%s' does not appear to be a muon bundle energy distribution", bundles.c_str());
}

double
SplineEnergyDistribution::operator()(double depth, double cos_theta, 
    unsigned multiplicity, double radius, double energy) const
{
	double coords[5] = {cos_theta, depth, multiplicity, radius, std::log(energy)};
	double logprob;
	
	if (multiplicity < 2) {
		coords[2] = coords[4];
		if (singles_.Eval(coords, &logprob) != 0)
			return 0.;
	} else if (bundles_.Eval(coords, &logprob) != 0)
		return 0.;
	
	return std::exp(logprob);
}

Sample
SplineEnergyDistribution::Generate(double depth, double cos_theta, 
    unsigned multiplicity, double radius) const
{
	log_fatal("Sampling is not yet implemented");
	return Sample(1., 1.);
}

OffsetPowerLaw::OffsetPowerLaw(double gamma, double offset, double emin, double emax)
    : gamma_(gamma), offset_(offset), emin_(emin), emax_(emax)
{
	nmin_ = std::pow(emin + offset, 1-gamma);
	nmax_ = std::pow(emax + offset, 1-gamma);
	norm_ = (1-gamma)/(nmax_ - nmin_);
}

double
OffsetPowerLaw::operator()(double energy) const
{
	return norm_*std::pow(energy + offset_, -gamma_);
}

double
OffsetPowerLaw::Generate() const
{
	return std::pow(rng_->Uniform()*(nmax_ - nmin_) + nmin_, 1./(1.-gamma_)) - offset_;
}

}