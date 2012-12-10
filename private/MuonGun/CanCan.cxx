
#include <MuonGun/CanCan.h>
#include <MuonGun/I3MuonGun.h>
#include <dataclasses/I3Constants.h>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <icetray/I3Module.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

namespace I3MuonGun {

void
BundleInjector::SetRandomService(I3RandomServicePtr r)
{
	Distribution::SetRandomService(r);
	flux_->SetRandomService(r);
	multiplicity_->SetRandomService(r);
}

void
BundleInjector::SetSurface(CylinderPtr p)
{
	surface_ = p;
	totalRate_ = NAN;
	CalculateMaxFlux();
}

void
BundleInjector::SetFlux(SingleMuonFluxPtr p)
{
	flux_ = p;
	totalRate_ = NAN;
	CalculateMaxFlux();
}

void
BundleInjector::SetMultiplicity(MultiplicityFractionPtr p)
{
	multiplicity_ = p;
	totalRate_ = NAN;
}

void
BundleInjector::CalculateMaxFlux()
{
	if (surface_ && flux_)
		maxFlux_ = (*flux_)(surface_->GetMinDepth(), 1.)*surface_->GetMaxDifferentialArea();
}

double
BundleInjector::GetTotalRate() const
{
	if (std::isnan(totalRate_) && surface_ && flux_ && multiplicity_) {
		totalRate_ = 0;
		for (unsigned m = multiplicity_->GetMin(); m <= multiplicity_->GetMax(); m++)
			totalRate_ += surface_->IntegrateFlux(detail::multiply<2>(
			    boost::cref(*flux_), boost::bind(boost::cref(*multiplicity_), _1, _2, m)));
	}
	return totalRate_;
}

boost::tuple<I3Particle, unsigned, double>
BundleInjector::GenerateAxis() const
{
	I3Direction dir;
	I3Position pos;
	unsigned m;
	double h, flux;
	do {
		surface_->SampleImpactRay(pos, dir, *rng_);
		m = rng_->Integer(multiplicity_->GetMax()-multiplicity_->GetMin())
		    + multiplicity_->GetMin();
		// Now, calculate the flux expectation at the chosen zenith angle
		// and at the depth where the shower axis crosses the surface
		double h = GetDepth(pos.GetZ());
		double coszen = cos(dir.GetZenith());
		flux = (*flux_)(h, coszen)
		    * (*multiplicity_)(h, coszen, m)
		    * surface_->GetDifferentialArea(coszen);
	} while (flux <= rng_->Uniform(0., maxFlux_));
	
	I3Particle p;
	p.SetPos(pos);
	p.SetDir(dir);
	p.SetShape(I3Particle::Primary);
	p.SetLocationType(I3Particle::Anywhere);
	p.SetType(I3Particle::unknown);
	
	return ::boost::make_tuple(p, m, GetTotalRate());
}


class BundleGenerator : public I3Module, private BundleInjector {
public:
	BundleGenerator(const I3Context &ctx) : I3Module(ctx), nEvents_(0)
	{
		AddOutBox("OutBox");
		AddParameter("NEvents", "", 1u);
	}
	
	void Configure()
	{
		BundleInjector::SetFlux(boost::make_shared<AdHocSingleMuonFlux>());
		MultiplicityFractionPtr multiplicity = boost::make_shared<BMSSMultiplicityFraction>();
		multiplicity->SetMin(1);
		multiplicity->SetMax(1);
		BundleInjector::SetMultiplicity(multiplicity);
		BundleInjector::SetSurface(boost::make_shared<Cylinder>(1600, 800));
		
		energyGenerator_ = boost::make_shared<OffsetPowerLaw>(2, 500., 50, 1e6);
		
		std::string singles, bundles;
		{
			std::ostringstream tablePath;
			tablePath << getenv("I3_BUILD") << "/MuonGun/resources/scripts/fitting/Hoerandel5_atmod12_SIBYLL.fits";
			singles = tablePath.str();
		}
		{
			std::ostringstream tablePath;
			tablePath << getenv("I3_BUILD") << "/MuonGun/resources/scripts/fitting/Hoerandel5_atmod12_SIBYLL_multi.fits";
			bundles = tablePath.str();
		}
		energySpectrum_ = boost::make_shared<SplineEnergyDistribution>(singles, bundles);
		
		std::string radial;
		{
			std::ostringstream tablePath;
			tablePath << getenv("I3_BUILD") << "/MuonGun/resources/scripts/fitting/Hoerandel5_atmod12_SIBYLL.radius.fits";
			radial = tablePath.str();
		}
		radialDistribution_ = boost::make_shared<SplineRadialDistribution>(radial);
		
		GetParameter("NEvents", maxEvents_);
		I3RandomServicePtr rng = context_.Get<I3RandomServicePtr>();
		BundleInjector::SetRandomService(rng);
		energyGenerator_->SetRandomService(rng);
	}
	
	void DAQ(I3FramePtr frame)
	{
		boost::tuple<I3Particle, unsigned, double> axis = BundleInjector::GenerateAxis();
		
		I3Particle &primary = boost::get<0>(axis);
		I3MCTreePtr mctree = boost::make_shared<I3MCTree>();
		I3MCTreeUtils::AddPrimary(mctree, primary);
		
		double h = GetDepth(primary.GetPos().GetZ());
		double coszen = cos(primary.GetDir().GetZenith());
		double rate = boost::get<2>(axis)/maxEvents_;
		
		unsigned multiplicity = boost::get<1>(axis);
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
			rate *= (*energySpectrum_)(h, coszen, multiplicity, radius, track.GetEnergy())
			    / (*energyGenerator_)(track.GetEnergy());
			I3MCTreeUtils::AppendChild(mctree, primary, track);
		}
		
		frame->Put("I3MCTree", mctree);
		frame->Put("Weight", boost::make_shared<I3Double>(rate));
		
		PushFrame(frame);
		
		if (++nEvents_ >= maxEvents_)
			RequestSuspension();
	}
private:
	size_t nEvents_, maxEvents_;
	boost::shared_ptr<OffsetPowerLaw> energyGenerator_;
	EnergyDistributionPtr energySpectrum_;
	RadialDistributionPtr radialDistribution_;

};

}

I3_MODULE(I3MuonGun::BundleGenerator);
