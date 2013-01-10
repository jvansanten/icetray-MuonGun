/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>
#include <icetray/I3Units.h>

namespace I3MuonGun {

RadialDistribution::~RadialDistribution() {}

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
BMSSRadialDistribution::operator()(double depth, double cos_theta,
    unsigned N, double radius) const
{
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double theta = acos(cos_theta);
	
	return GetGenerationProbability(GetMeanRadius(h, theta, N), GetShapeParameter(h, theta, N), radius);
}

double
BMSSRadialDistribution::Generate(I3RandomService &rng, double depth, double cos_theta,
    unsigned N) const
{
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
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
		r = rng.Uniform(rmax_);
	} while (rng.Uniform(max_prob) <= GetGenerationProbability(R, a, r));

	return r;
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

double
SplineRadialDistribution::Generate(I3RandomService &rng, double depth,
    double cos_theta, unsigned N) const
{
	double radius, logprob, maxprob;
	std::pair<double, double> extent = I3SplineTable::GetExtents(3);
	double coords[4] = {cos_theta, depth, N, extent.first};
	if (I3SplineTable::Eval(coords, &maxprob) != 0)
		maxprob = -std::numeric_limits<double>::infinity();
	
	// The spline is fit to log(dP/dr^2) as a function of r,
	// so we generate proposals uniformly in r^2, then take
	// a square root to evaluate.
	do {
		coords[3] = std::sqrt(rng.Uniform(extent.first*extent.first,
		    extent.second*extent.second));
		if (I3SplineTable::Eval(coords, &logprob) != 0)
			logprob = -std::numeric_limits<double>::infinity();
	} while (std::log(rng.Uniform()) > logprob - maxprob);
	
	return coords[3];
}

}
