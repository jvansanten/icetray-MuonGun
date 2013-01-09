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
	tablePath << getenv("I3_BUILD") << "/MuonGun/resources/scripts/fitting/" << subpath;
	return tablePath.str();
}

}

EnergyDependentSurfaceInjector::EnergyDependentSurfaceInjector()
{
	surface_ = boost::make_shared<Cylinder>(1600, 800);
	
	flux_ = boost::make_shared<SplineFlux>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_flux.fits"),
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_flux.fits"));
	flux_->SetMinMultiplicity(1);
	flux_->SetMaxMultiplicity(1);
	
	energyGenerator_ = boost::make_shared<OffsetPowerLaw>(2, 500., 50, 1e6);
	
	radialDistribution_ = boost::make_shared<SplineRadialDistribution>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.radius.fits"));
}

GenerationProbabilityPtr
EnergyDependentSurfaceInjector::Clone() const
{
	return boost::make_shared<EnergyDependentSurfaceInjector>(*this);
}

void
EnergyDependentSurfaceInjector::Generate(I3RandomService &rng, I3MCTree &tree,
    BundleConfiguration &bundle) const
{
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double flux, max_flux;
	double max_energy;
	double h, coszen;
	SamplingSurfaceConstPtr surface;
	do {
		// Choose a multiplicity
		bundle.clear();
		m = rng.Integer(flux_->GetMaxMultiplicity() - flux_->GetMinMultiplicity())
		    + flux_->GetMinMultiplicity();
		// Choose an ensemble of energies
		max_energy = 0.;
		for (unsigned i=0; i < m; i++) {
			bundle.push_back(std::make_pair(0., energyGenerator_->Generate(rng)));
			if (bundle.back().second > max_energy)
				max_energy = bundle.back().second;
		}
		// Choose sampling surface based on highest-energy muon
		surface = GetSurface(max_energy);
		// Sample from energy-dependent surface
		surface->SampleImpactRay(pos, dir, rng);
		// Calculate the differential flux expectation for
		// this surface at the chosen angle, depth, and
		// multiplicity, and compare to the maximum differential
		// flux anywhere on the surface.
		h = GetDepth(pos.GetZ());
		coszen = cos(dir.GetZenith());
		flux = (*flux_)(h, coszen, m)
		    * surface->GetDifferentialArea(coszen);
		max_flux = (*flux_)(surface->GetMinDepth(),
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
		track.SetEnergy(bspec.second);
		bspec.first = radius;
		I3MCTreeUtils::AppendChild(tree, primary, track);
	}
}

SamplingSurfacePtr
EnergyDependentSurfaceInjector::GetSurface(double energy) const
{
	// TODO implement realistic/configurable scaling
	CylinderPtr surface = boost::make_shared<Cylinder>(*surface_);
	
	I3Position center = surface->GetCenter();
	double ebounds[2] = {1e3, 1e4};
	double rbounds[2] = {125, 800};
	double tbounds[2] = {-100, 800};
	if (energy < ebounds[0])
		surface->SetRadius(rbounds[0]);
	else if (energy > ebounds[1])
		surface->SetRadius(rbounds[1]);
	else
		surface->SetRadius(rbounds[0] + (rbounds[1]-rbounds[0])*(energy - ebounds[0])/(ebounds[1]-ebounds[0]));
	
	double top;
	if (energy < ebounds[0])
		top = tbounds[0];
	else if (energy > ebounds[1])
		top = tbounds[1];
	else
		top = tbounds[0] + (tbounds[1]-tbounds[0])*(energy - ebounds[0])/(ebounds[1]-ebounds[0]);
	surface->SetLength(top+800.);
	center.SetZ(top-surface->GetLength()/2);
	surface->SetCenter(center);
	
	return surface;
}

SamplingSurfaceConstPtr
EnergyDependentSurfaceInjector::GetInjectionSurface(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	double energy = 0.;
	BOOST_FOREACH(const BundleConfiguration::value_type &bspec, bundle)
		if (bspec.second > energy)
			energy = bspec.second;
	
	return GetSurface(energy);
}

double
EnergyDependentSurfaceInjector::GetTotalRate(double energy) const
{
	SamplingSurfaceConstPtr surface = GetSurface(energy);
	
	return surface->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, 1u));
}

double
EnergyDependentSurfaceInjector::GetGenerationProbability(double h,
    double coszen, const BundleConfiguration &bundle) const
{
	unsigned m = bundle.size();
	double prob = flux_->operator()(h, coszen, m);
	double max_energy = 0.;
	BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundle) {
		if (m > 1)
			prob *= (*radialDistribution_)(h, coszen, m, pair.first);
		prob *= (*energyGenerator_)(pair.second);
		if (pair.second > max_energy)
			max_energy = pair.second;
	}
	// FIXME: GetTotalRate() is potentially expensive, as are repeated heap
	// allocations in GetSurface()
	prob *= GetSurface(max_energy)->GetDifferentialArea(coszen)/GetTotalRate(max_energy);
	
	return prob;
}

}
