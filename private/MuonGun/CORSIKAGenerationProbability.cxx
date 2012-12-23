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
CORSIKAGenerationProbability::GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// assert(steps.first >= 0);
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundlespec.size();
	
	double prob = flux_->operator()(h, coszen, m)*surface_->GetDifferentialArea(coszen);
	BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundlespec) {
		if (m > 1)
			prob *= (*radialDistribution_)(h, coszen, m, pair.first);
		prob *= (*energyDistribution_)(h, coszen, m, pair.first, pair.second);
	}
	
	return prob;
}

}
