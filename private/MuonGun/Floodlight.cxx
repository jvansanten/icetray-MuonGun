
#include <MuonGun/Floodlight.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Surface.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

#include <boost/foreach.hpp>

namespace I3MuonGun {

Floodlight::Floodlight() : surface_(new Cylinder(1000, 600, I3Position(31.25, 19.64, 0))),
    energyGenerator_(new OffsetPowerLaw(1.1, 0., 5e2, 1e7)), totalRate_(NAN)
{}

GenerationProbabilityPtr
Floodlight::Clone() const
{
	return GenerationProbabilityPtr(new Floodlight(*this));
}

void
Floodlight::Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const
{
	I3Direction dir;
	I3Position pos;
	surface_->SampleImpactRay(pos, dir, rng, -1, 1);
	
	I3Particle primary;
	primary.SetDir(dir);
	primary.SetPos(pos);
	primary.SetTime(0);
	primary.SetEnergy(energyGenerator_->Generate(rng));
	
	I3Particle muon(primary);
	muon.SetType(I3Particle::MuMinus);
	muon.SetLocationType(I3Particle::InIce);
	muon.SetShape(I3Particle::Null);
	
	I3MCTreeUtils::AddPrimary(tree, primary);
	I3MCTreeUtils::AppendChild(tree, primary, muon);
}

static double
one(double,double) { return 1.; }

double
Floodlight::GetTotalRate() const
{
	if (std::isfinite(totalRate_))
		return totalRate_;
	totalRate_ = std::log(surface_->IntegrateFlux(one, -1, 1));
	return totalRate_;
}

double
Floodlight::GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// This shower axis doesn't intersect the sampling surface. Bail.
	if (!std::isfinite(steps.first))
		return -std::numeric_limits<double>::infinity();
	
	double h = GetDepth(axis.GetPos().GetZ() + steps.first*axis.GetDir().GetZ());
	double coszen = cos(axis.GetDir().GetZenith());
	unsigned m = bundle.size();
	double logprob = std::log(surface_->GetDifferentialArea(coszen));
	BOOST_FOREACH(const BundleEntry &track, bundle)
		logprob += energyGenerator_->GetLog(track.energy);
	
	return logprob - GetTotalRate();
}

}