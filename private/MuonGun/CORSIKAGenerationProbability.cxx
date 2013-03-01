/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/I3MuonGun.h>
#include <MuonGun/CORSIKAGenerationProbability.h>

#include <dataclasses/physics/I3Particle.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <boost/foreach.hpp>

namespace I3MuonGun {

CORSIKAGenerationProbability::CORSIKAGenerationProbability(SamplingSurfacePtr s, FluxPtr f, RadialDistributionPtr r, EnergyDistributionPtr e)
    : surface_(s), flux_(f), radialDistribution_(r), energyDistribution_(e)
{
	if (!surface_)
		log_fatal("No sampling surface defined!");
	if (!flux_)
		log_fatal("No flux defined!");
	if (!radialDistribution_)
		log_fatal("No radial distribution defined!");
	if (!energyDistribution_)
		log_fatal("No energy distribution defined!");
}

GenerationProbabilityPtr
CORSIKAGenerationProbability::Clone() const
{
	return boost::make_shared<CORSIKAGenerationProbability>(*this);
}

double
CORSIKAGenerationProbability::GetLogBundleGenerationProbability(const I3Particle &axis,
    unsigned m) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return -std::numeric_limits<double>::infinity();
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	return flux_->GetLog(h, coszen, m) + std::log(surface_->GetDifferentialArea(coszen));
}

double
CORSIKAGenerationProbability::GetLogGenerationProbability(double h, double coszen,
    unsigned m, double radius, double energy) const
{
	return energyDistribution_->GetLog(h, coszen, m, radius, energy);
}

double
CORSIKAGenerationProbability::GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundlespec) const
{	
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return 0.;
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundlespec.size();
	double logprob = flux_->GetLog(h, coszen, m) + std::log(surface_->GetDifferentialArea(coszen));
	BOOST_FOREACH(const BundleConfiguration::value_type &track, bundlespec) {
		if (m > 1)
			logprob += radialDistribution_->GetLog(h, coszen, m, track.radius);
		logprob += energyDistribution_->GetLog(h, coszen, m, track.radius, track.energy);
	}
	
	return logprob;
}

}
