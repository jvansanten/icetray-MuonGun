
#include <MuonGun/CanCan.h>
#include <MuonGun/I3MuonGun.h>
#include <dataclasses/I3Constants.h>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

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

StaticSurfaceInjector::StaticSurfaceInjector()
{
	SetSurface(boost::make_shared<Cylinder>(1600, 800));
	
	FluxPtr flux = boost::make_shared<SplineFlux>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.single_flux.fits"),
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_flux.fits"));
	flux->SetMinMultiplicity(1);
	flux->SetMaxMultiplicity(1);
	SetFlux(flux);
	
	energyGenerator_ = boost::make_shared<OffsetPowerLaw>(2, 500., 50, 1e6);
	
	// energySpectrum_ = boost::make_shared<SplineEnergyDistribution>(
	//     GetTablePath("Hoerandel5_atmod12_SIBYLL.single_energy.fits"),
	//     GetTablePath("Hoerandel5_atmod12_SIBYLL.bundle_energy.fits"));
	
	radialDistribution_ = boost::make_shared<SplineRadialDistribution>(
	    GetTablePath("Hoerandel5_atmod12_SIBYLL.radius.fits"));
}

GenerationProbabilityPtr
StaticSurfaceInjector::Clone() const
{
	return boost::make_shared<StaticSurfaceInjector>(*this);
}

void
StaticSurfaceInjector::SetRandomService(I3RandomServicePtr r)
{
	Distribution::SetRandomService(r);
	flux_->SetRandomService(r);
	energyGenerator_->SetRandomService(r);
	radialDistribution_->SetRandomService(r);
}

void
StaticSurfaceInjector::SetSurface(SamplingSurfacePtr p)
{
	surface_ = p;
	totalRate_ = NAN;
	CalculateMaxFlux();
}

void
StaticSurfaceInjector::SetFlux(FluxPtr p)
{
	flux_ = p;
	totalRate_ = NAN;
	CalculateMaxFlux();
}

void
StaticSurfaceInjector::CalculateMaxFlux()
{
	if (surface_ && flux_)
		maxFlux_ = (*flux_)(surface_->GetMinDepth(), 1., 1u)*surface_->GetMaxDifferentialArea();
}

double
StaticSurfaceInjector::GetTotalRate() const
{
	if (std::isnan(totalRate_) && surface_ && flux_) {
		totalRate_ = 0;
		for (unsigned m = flux_->GetMinMultiplicity(); m <= flux_->GetMaxMultiplicity(); m++) {
			std::cout << "m = " << m << std::endl;
			totalRate_ += surface_->IntegrateFlux(boost::bind(boost::cref(*flux_), _1, _2, m));
		}
	}
	return totalRate_;
}

void
StaticSurfaceInjector::GenerateAxis(std::pair<I3Particle, unsigned> &axis) const
{
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double flux;
	do {
		surface_->SampleImpactRay(pos, dir, *rng_);
		m = rng_->Integer(flux_->GetMaxMultiplicity() - flux_->GetMinMultiplicity())
		    + flux_->GetMinMultiplicity();
		// Now, calculate the flux expectation at the chosen zenith angle
		// and at the depth where the shower axis crosses the surface
		double h = GetDepth(pos.GetZ());
		double coszen = cos(dir.GetZenith());
		flux = (*flux_)(h, coszen, m)
		    * surface_->GetDifferentialArea(coszen);
	} while (flux <= rng_->Uniform(0., maxFlux_));
	
	axis.first.SetPos(pos);
	axis.first.SetDir(dir);
	axis.first.SetShape(I3Particle::Primary);
	axis.first.SetLocationType(I3Particle::Anywhere);
	axis.first.SetType(I3Particle::unknown);
	axis.second = m;
}

void
StaticSurfaceInjector::FillMCTree(const std::pair<I3Particle, unsigned> &axis,
    I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	const I3Particle &primary = axis.first;
	I3MCTreeUtils::AddPrimary(mctree, primary);
	double h = GetDepth(primary.GetPos().GetZ());
	double coszen = cos(primary.GetDir().GetZenith());
	
	unsigned multiplicity = axis.second;
	for (unsigned i=0; i < multiplicity; i++) {
		I3Particle track;
		track.SetLocationType(I3Particle::InIce);
		track.SetType(I3Particle::MuMinus);
		track.SetDir(primary.GetDir());
		track.SetSpeed(I3Constants::c);
		
		double radius = (multiplicity > 1u) ?
		    radialDistribution_->Generate(h, coszen, multiplicity).value : 0.;
		double azimuth = rng_->Uniform(0, 2*M_PI);
		I3Position offset(radius, 0, 0);
		offset.RotateY(primary.GetDir().GetZenith());
		offset.RotateZ(azimuth);
		track.SetPos(offset.GetX()+primary.GetPos().GetX(),
		             offset.GetY()+primary.GetPos().GetY(),
		             offset.GetZ()+primary.GetPos().GetZ());
		// TODO: shift the times and positions of off-axis tracks
		// so that they correspond to a plane wave crossing the sampling
		// surface at time 0
			
		track.SetEnergy(energyGenerator_->Generate());
		I3MCTreeUtils::AppendChild(mctree, primary, track);
		bundlespec.push_back(std::make_pair(radius, track.GetEnergy()));
	}
}

void
StaticSurfaceInjector::Generate(I3MCTree &mctree, BundleConfiguration &bundlespec) const
{
	std::pair<I3Particle, unsigned> axis;
	GenerateAxis(axis);
	FillMCTree(axis, mctree, bundlespec);
}

double
StaticSurfaceInjector::GetGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundlespec) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// assert(steps.first >= 0);
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundlespec.size();
	
	double prob = flux_->operator()(h, coszen, m)*surface_->GetDifferentialArea(coszen)/GetTotalRate();
	BOOST_FOREACH(const BundleConfiguration::value_type &pair, bundlespec) {
		if (m > 1)
			prob *= (*radialDistribution_)(h, coszen, m, pair.first);
		prob *= (*energyGenerator_)(pair.second);
	}
	
	return prob;
}

}

// I3_MODULE(I3MuonGun::BundleGenerator);
