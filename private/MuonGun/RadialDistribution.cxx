
#include <MuonGun/RadialDistribution.h>
#include <icetray/I3Units.h>

namespace I3MuonGun {

template <typename Archive>
void
RadialDistribution::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("Distribution", base_object<Distribution>(*this));
}

template <typename Archive>
void
BMSSRadialDistribution::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("RadialDistribution", base_object<RadialDistribution>(*this));
	ar & make_nvp("Rho0a",   rho0a_);
	ar & make_nvp("Rho0b",   rho0b_);
	ar & make_nvp("Rho1",    rho1_);
	ar & make_nvp("Theta0",  theta0_);
	ar & make_nvp("F",       f_);
	ar & make_nvp("Alpha0a", alpha0a_);
	ar & make_nvp("Alpha0b", alpha0b_);
	ar & make_nvp("Alpha1a", alpha1a_);
	ar & make_nvp("Alpha1b", alpha1b_);
	ar & make_nvp("RMax",    rmax_);
}

BMSSRadialDistribution::BMSSRadialDistribution() : rho0a_(-1.786), rho0b_(28.26),
    rho1_(-1.06), theta0_(1.3), f_(10.4), alpha0a_(-0.448), alpha0b_(4.969),
    alpha1a_(0.0194), alpha1b_(0.276), rmax_(250*I3Units::m) {};

double
BMSSRadialDistribution::GetMeanRadius(double h, double theta, unsigned N) const
{
	return ((rho0a_*N + rho0b_)*pow(h,rho1_))/(exp((theta-theta0_)*f_)+1.);
}

double
BMSSRadialDistribution::GetShapeParameter(double h, double theta, unsigned N) const
{
	return (alpha0a_*N + alpha0b_)*exp(h*(alpha1a_*N + alpha1b_));
}

double
BMSSRadialDistribution::GetGenerationProbability(double R, double a, double radius) const
{
	double R0 = R*(a-3)/2.;
	
	return (a-1)*(a-2)*pow(R0, a-2)*(radius/pow(radius+R0, a));
}

double
BMSSRadialDistribution::operator()(double h, double cos_theta,
    unsigned N, double radius) const
{
	double theta = acos(cos_theta);
	
	return GetGenerationProbability(GetMeanRadius(h, theta, N), GetShapeParameter(h, theta, N), radius);
}

Sample
BMSSRadialDistribution::Generate(double h, double cos_theta,
    unsigned N) const
{
	if (!rng_)
		log_fatal("No random number service set!");

	double theta = acos(cos_theta);
	double R = GetMeanRadius(h, theta, N);
	double a = GetShapeParameter(h, theta, N);
	
	double peak_radius = R*(a-3)/(2.*(a-1));
	if (!std::isfinite(peak_radius))
		log_fatal("Peak radius is not finite!");
	double max_prob = GetGenerationProbability(R, a, peak_radius);
	if (!std::isfinite(max_prob))
		log_fatal("Peak probability is not finite!");
	double r;
	do {
		r = rng_->Uniform(rmax_);
	} while (rng_->Uniform(max_prob) <= GetGenerationProbability(R, a, r));

	return Sample(r, GetGenerationProbability(R, a, r));
}

SplineRadialDistribution::SplineRadialDistribution(const std::string &path)
    : I3SplineTable(path) {}

double
SplineRadialDistribution::operator()(double depth, double cos_theta,
    unsigned N, double radius) const
{
	double coords[4] = {cos_theta, depth, N, radius};
	double logprob;
	
	if (I3SplineTable::Eval(coords, &logprob) != 0)
		return 0.;
	else
		// Spline is fit to log(dP/dr^2)
		return 2*radius*std::exp(logprob);
}

Sample
SplineRadialDistribution::Generate(double depth, double cos_theta,
    unsigned N) const
{
	if (!rng_)
		log_fatal("No random number service set!");
	
	double radius, logprob, maxprob;
	std::pair<double, double> extent = I3SplineTable::GetExtents(3);
	double coords[4] = {cos_theta, depth, N, extent.first};
	if (I3SplineTable::Eval(coords, &maxprob) != 0)
		maxprob = -std::numeric_limits<double>::infinity();
	
	// The spline is fit to log(dP/dr^2) as a function of r,
	// so we generate proposals uniformly in r^2, then take
	// a square root to evaluate.
	do {
		coords[3] = std::sqrt(rng_->Uniform(extent.first*extent.first,
		    extent.second*extent.second));
		if (I3SplineTable::Eval(coords, &logprob) != 0)
			logprob = -std::numeric_limits<double>::infinity();
	} while (std::log(rng_->Uniform()) > logprob - maxprob);
	
	return Sample(coords[3], 2*coords[3]*std::exp(logprob));
}

}

// I3_SERIALIZABLE(I3MuonGun::RadialDistribution);
// I3_SERIALIZABLE(I3MuonGun::BMSSRadialDistribution);


