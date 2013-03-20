/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Surface.h>
#include <MuonGun/Flux.h>
#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/RadialDistribution.h>
#include <dataclasses/I3Constants.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/I3RandomService.h>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

namespace I3MuonGun {

namespace {

std::string
GetTablePath(const std::string &subpath)
{
	std::ostringstream tablePath;
	tablePath << getenv("I3_BUILD") << "/MuonGun/resources/tables/" << subpath;
	return tablePath.str();
}

SamplingSurfacePtr
ScaleForIC79(double energy)
{
	// Distance to border of detector scales as (b*(a - log10(e)))**(1./alpha)
	double a = 4.5;
	double b = 1e5;
	int alpha = 2;
	
	// Center of cylinder moves from barycenter of total IC79 geometry at high energies
	// to the center of DeepCore (String 36) at low energies
	double x[2] = { 46.29, 31.25};
	double y[2] = {-34.88, 19.64};
	
	double d = log10(energy) >= a ? 0 : pow(b*(a - log10(energy)), 1./alpha);
	double r = std::max(600. - d, 100.);
	double l = std::max(1000. - d, 400.);
	I3Position center(x[0] + (r/600)*(x[1]-x[0]), y[0] + (r/600)*(y[1]-y[0]), -500 + l/2.);
	
	return boost::make_shared<Cylinder>(l, r, center);
}

}

EnergyDependentSurfaceInjector::EnergyDependentSurfaceInjector(FluxPtr flux, RadialDistributionPtr radius,
    boost::shared_ptr<OffsetPowerLaw> energies, boost::function<SamplingSurfacePtr (double)> scaling)
    : scalingFunction_(scaling), flux_(flux), energyGenerator_(energies), radialDistribution_(radius)
{
	if (!flux_) {
		flux_ = boost::make_shared<SplineFlux>(
		    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_flux.fits"),
		    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_flux.fits"));
		flux_->SetMinMultiplicity(1);
		flux_->SetMaxMultiplicity(1);
	}
	if (!radialDistribution_)
		radialDistribution_ = boost::make_shared<SplineRadialDistribution>(
		    GetTablePath("Hoerandel5_atmod12_SIBYLL.radius.fits"));
	if (!energyGenerator_)
		energyGenerator_ = boost::make_shared<OffsetPowerLaw>(2, 500., 50, 1e6);
	if (!scalingFunction_)
		scalingFunction_ = &ScaleForIC79;
	
	injectionSurface_ = boost::make_shared<Cylinder>(1600, 800);
}

GenerationProbabilityPtr
EnergyDependentSurfaceInjector::Clone() const
{
	return boost::make_shared<EnergyDependentSurfaceInjector>(*this);
}

bool
EnergyDependentSurfaceInjector::IsCompatible(GenerationProbabilityConstPtr o) const
{
	boost::shared_ptr<const EnergyDependentSurfaceInjector> other
	    = boost::dynamic_pointer_cast<const EnergyDependentSurfaceInjector>(o);
	if (!other)
		return false;
	else
		return (*flux_ == *(other->flux_)
		    && *radialDistribution_ == *(other->radialDistribution_)
		    && *energyGenerator_ == *(other->energyGenerator_)
		    && (scalingFunction_ == &ScaleForIC79) && (other->scalingFunction_ == &ScaleForIC79));
}

void
EnergyDependentSurfaceInjector::Generate(I3RandomService &rng, I3MCTree &tree,
    BundleConfiguration &bundle) const
{
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double flux, max_flux;
	double h, coszen;
	SamplingSurfaceConstPtr surface;
	do {
		// Choose a multiplicity
		bundle.clear();
		m = rng.Integer(flux_->GetMaxMultiplicity() - flux_->GetMinMultiplicity())
		    + flux_->GetMinMultiplicity();
		// Choose an ensemble of energies
		for (unsigned i=0; i < m; i++)
			bundle.push_back(BundleEntry(0., energyGenerator_->Generate(rng)));
		bundle.sort();
		// Choose target surface based on highest-energy muon
		surface = GetTargetSurface(bundle.front().energy);
		// Sample an impact point on the target surface
		surface->SampleImpactRay(pos, dir, rng);
		// Snap the impact point back to the injection surface
		std::pair<double, double> steps = injectionSurface_->GetIntersection(pos, dir);
		if (!(steps.first <= 0))
			log_fatal("The target point is outside the injection surface!");
		pos.SetX(pos.GetX() + steps.first*dir.GetX());
		pos.SetY(pos.GetY() + steps.first*dir.GetY());
		pos.SetZ(pos.GetZ() + steps.first*dir.GetZ());
		// Calculate the differential flux expectation for
		// this surface at the chosen angle, depth, and
		// multiplicity, and compare to the maximum differential
		// flux anywhere on the surface.
		h = GetDepth(pos.GetZ());
		coszen = cos(dir.GetZenith());
		flux = (*flux_)(h, coszen, m)
		    * surface->GetDifferentialArea(coszen);
		max_flux = (*flux_)(injectionSurface_->GetMinDepth(),
		    1., 1u)*surface->GetMaxDifferentialArea();
	} while (flux <= rng.Uniform(0., max_flux));
	
	I3Particle primary;
	primary.SetPos(pos);
	primary.SetDir(dir);
	primary.SetShape(I3Particle::Primary);
	primary.SetLocationType(I3Particle::Anywhere);
	primary.SetType(I3Particle::unknown);
	primary.SetTime(0.);
	I3MCTreeUtils::AddPrimary(tree, primary);
	
	// For each muon, draw a radial offset and add an entry
	// to the MCTree
	BOOST_FOREACH(BundleConfiguration::value_type &bspec, bundle) {
		double radius = 0., azimuth = 0.;
		if (m > 1u) {
			radius = radialDistribution_->Generate(rng, h, coszen, m);
			azimuth = rng.Uniform(0., 2*M_PI);
		}
		
		I3Particle track = CreateParallelTrack(radius, azimuth, *surface, primary);
		track.SetEnergy(bspec.energy);
		bspec.radius = radius;
		I3MCTreeUtils::AppendChild(tree, primary, track);
	}
}

SamplingSurfacePtr
EnergyDependentSurfaceInjector::GetTargetSurface(double energy) const
{
	if (scalingFunction_)
		return scalingFunction_(energy);
	else
		return boost::make_shared<Cylinder>(1600, 800);
}

SamplingSurfaceConstPtr
EnergyDependentSurfaceInjector::GetInjectionSurface(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	return injectionSurface_;
}

double
EnergyDependentSurfaceInjector::GetTotalRate(SamplingSurfaceConstPtr surface) const
{
	double total_rate = 0.;
	for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++)
		total_rate += surface->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, m));
	return total_rate;
}

double
EnergyDependentSurfaceInjector::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	// Entries are sorted in descending order of energy, so the
	// "minimum" entry has the maximum energy
	SamplingSurfaceConstPtr surface =
	    GetTargetSurface(std::min_element(bundle.begin(), bundle.end())->energy);
	std::pair<double, double> steps =
	    surface->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return 0.;
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundle.size();
	double logprob = flux_->GetLog(h, coszen, m);
	double max_energy = 0.;
	BOOST_FOREACH(const BundleEntry &track, bundle) {
		if (m > 1)
			logprob += radialDistribution_->GetLog(h, coszen, m, track.radius);
		logprob += energyGenerator_->GetLog(track.energy);
	}
	// FIXME: rate integration is potentially expensive, as are repeated heap
	// allocations in GetSurface()
	logprob += std::log(surface->GetDifferentialArea(coszen)) - std::log(GetTotalRate(surface));
	
	return logprob;
}

}
