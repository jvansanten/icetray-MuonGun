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
	BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundlespec) {
		if (m > 1)
			logprob += radialDistribution_->GetLog(h, coszen, m, pair.first);
		logprob += energyDistribution_->GetLog(h, coszen, m, pair.first, pair.second);
	}
	
	return logprob;
}

}
