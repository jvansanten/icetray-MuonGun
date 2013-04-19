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

}

EnergyDependentSurfaceInjector::EnergyDependentSurfaceInjector(FluxPtr flux, RadialDistributionPtr radius,
    boost::shared_ptr<OffsetPowerLaw> energies, SurfaceScalingFunctionPtr scaling)
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
		    && *scalingFunction_ == *(other->scalingFunction_));
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
		return scalingFunction_->GetSurface(energy);
	else
		return boost::make_shared<Cylinder>(1600, 800);
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
		return -std::numeric_limits<double>::infinity();
	
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
	// Bundle axes are sampled uniformly in the projected area of the target surface.
	// Here we apply a density correction to account for the fact that a locally isotropic
	// flux through the inner target surface is not necessarily an isotropic flux through
	// the outer surface when they have different shapes.
	double aspect_ratio = (injectionSurface_->GetDifferentialArea(coszen)/injectionSurface_->GetTotalArea())
	    / (surface->GetDifferentialArea(coszen)/surface->GetTotalArea());
	logprob += std::log(aspect_ratio) - std::log(GetTotalRate(surface));
	
	return logprob;
}

SurfaceScalingFunction::~SurfaceScalingFunction() {}

template <typename Archive>
void
SurfaceScalingFunction::serialize(Archive &ar, unsigned)
{}

template <typename Archive>
void
BasicSurfaceScalingFunction::serialize(Archive &ar, unsigned)
{
	ar & make_nvp("SurfaceScalingFunction", base_object<SurfaceScalingFunction>(*this));
}

BasicSurfaceScalingFunction::BasicSurfaceScalingFunction() :
    scale_(800., 240957.5), energyScale_(4., 4.), offset_(3.778, 3.622), power_(1.10, 2.23),
    rBounds_(0., 525.), zBounds_(-500., 400.),
    centerBounds_(pair(46.29,-34.88), pair(31.25, 19.64))
{}

BasicSurfaceScalingFunction::~BasicSurfaceScalingFunction() {}

double
BasicSurfaceScalingFunction::GetMargin(double logenergy, double scale, double offset, double power) const
{
	if (logenergy < offset)
		return pow(scale*(offset - logenergy), 1./power);
	else
		return 0.;
}

SamplingSurfacePtr
BasicSurfaceScalingFunction::GetSurface(double energy) const
{
	// Shrink the cylinder down by an energy-dependent amount
	double z = std::max(zBounds_.second -
	    GetMargin(std::log10(energy/energyScale_.first), scale_.first, offset_.first, power_.first), zBounds_.first);
	
	// Shrink the sides of the cylinder
	double r = std::max(rBounds_.second -
	    GetMargin(std::log10(energy/energyScale_.second), scale_.second, offset_.second, power_.second), rBounds_.first);
	
	// Move the center smoothly between the configured bounds
	double hscale = r/(rBounds_.second-rBounds_.first);
	I3Position center(centerBounds_.first.first + hscale*(centerBounds_.second.first-centerBounds_.first.first),
	    centerBounds_.first.second + hscale*(centerBounds_.second.second-centerBounds_.first.second),
	    (zBounds_.first + z)/2.);
	
	return boost::make_shared<Cylinder>(z-zBounds_.first, r, center);
}

void
BasicSurfaceScalingFunction::SetCapScaling(double energyScale, double scale, double offset, double power)
{
	energyScale_.first = energyScale;
	scale_.first = scale;
	offset_.first = offset;
	power_.first = power;
}

void
BasicSurfaceScalingFunction::SetSideScaling(double energyScale, double scale, double offset, double power)
{
	energyScale_.second = energyScale;
	scale_.second = scale;
	offset_.second = offset;
	power_.second = power;
}

void
BasicSurfaceScalingFunction::SetRadiusBounds(double rmin, double rmax)
{
	rBounds_ = std::make_pair(rmin, rmax);
}

void
BasicSurfaceScalingFunction::SetZBounds(double zmin, double zmax)
{
	zBounds_ = std::make_pair(zmin, zmax);
}

bool
BasicSurfaceScalingFunction::operator==(const SurfaceScalingFunction &o) const
{
	const BasicSurfaceScalingFunction *other = dynamic_cast<const BasicSurfaceScalingFunction*>(&o);
	if (!other)
		return false;
	else
		return (scale_ == other->scale_ &&
		    energyScale_ == other->energyScale_ &&
		    offset_ == other->offset_ &&
		    power_ == other->power_ &&
		    rBounds_ == other->rBounds_ &&
		    zBounds_ == other->zBounds_ &&
		    centerBounds_ == other->centerBounds_);
}

}
