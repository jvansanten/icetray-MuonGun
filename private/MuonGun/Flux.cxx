/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Flux.h>
#include <icetray/I3Units.h>
#include <limits>

namespace I3MuonGun {

Flux::Flux() : minMultiplicity_(1), maxMultiplicity_(1) {};

Flux::~Flux() {};

double
Flux::operator()(double depth, double cos_theta, unsigned multiplicity) const
{
	return std::exp(GetLog(depth, cos_theta, multiplicity));
}

BMSSFlux::BMSSFlux() : k0a_(7.2e-3), k0b_(-1.927), k1a_(-0.581), k1b_(0.034),
    v0a_(0.01041), v0b_(0.09912), v0c_(2.712), v1a_(0.01615), v1b_(0.6010)
{}
	
double
BMSSFlux::GetLog(double depth, double cos_theta, unsigned multiplicity) const
{
	// Convert to water-equivalent depth
	double h = (200*I3Units::m/I3Units::km)*0.832 + (depth-(200*I3Units::m/I3Units::km))*0.917;
	double flux = k0a_*std::pow(h, k0b_)*cos_theta*std::exp((k1a_*h + k1b_)/cos_theta);
	if (multiplicity > 1)
		flux *= std::pow(multiplicity,
		    -(v0a_*h*h + v0b_*h + v0c_)*std::exp(v1a_*std::exp(v1b_*h)/cos_theta));
	
	// Not the most stable thing, but only for demo purposes
	return std::log(flux);
}

SplineFlux::SplineFlux(const std::string &singles, const std::string &bundles)
    : singles_(singles), bundles_(bundles)
{
	SetMinMultiplicity(1);
	SetMaxMultiplicity(bundles_.GetExtents(2).second);
}

double
SplineFlux::GetLog(double depth, double cos_theta, unsigned multiplicity) const
{
	double coords[3] = {cos_theta, depth, multiplicity};
	double logflux;
	
	if (multiplicity < GetMinMultiplicity() || multiplicity > GetMaxMultiplicity())
		return -std::numeric_limits<double>::infinity();
	else if ((multiplicity > 1 ? bundles_ : singles_).Eval(coords, &logflux) != 0)
		return -std::numeric_limits<double>::infinity();
	else
		return logflux;
}


}