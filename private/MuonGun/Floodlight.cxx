
#include <MuonGun/Floodlight.h>
#include <MuonGun/I3MuonGun.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <MuonGun/Cylinder.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

#include <boost/foreach.hpp>

namespace I3MuonGun {

template <typename Archive>
void
Floodlight::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Generator", base_object<Generator>(*this));
	ar & make_nvp("Surface", surface_);
	ar & make_nvp("EnergySpectrum", energyGenerator_);
}

Floodlight::Floodlight(SamplingSurfacePtr surface, boost::shared_ptr<OffsetPowerLaw> energyGenerator) : surface_(surface),
    energyGenerator_(energyGenerator)
{
	if (!surface_)
		surface_ = boost::make_shared<Cylinder>(1000, 600, I3Position(31.25, 19.64, 0));
	if (!energyGenerator_)
		energyGenerator_ = boost::make_shared<OffsetPowerLaw>(1, 0., 5e2, 1e7);
}

GenerationProbabilityPtr
Floodlight::Clone() const
{
	return GenerationProbabilityPtr(new Floodlight(*this));
}

bool
Floodlight::IsCompatible(GenerationProbabilityConstPtr o) const
{
	boost::shared_ptr<const Floodlight> other = boost::dynamic_pointer_cast<const Floodlight>(o);
	if (!other)
		return false;
	else
		return (*surface_ == *(other->surface_)
		    && *energyGenerator_ == *(other->energyGenerator_));
}

void
Floodlight::Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle __attribute__((unused))) const
{
	I3Direction dir;
	I3Position pos;
	surface_->SampleImpactRay(pos, dir, rng, -1, 1);
	
	I3Particle primary;
	primary.SetDir(dir);
	primary.SetPos(pos);
	primary.SetTime(0);
	primary.SetEnergy(energyGenerator_->Generate(rng));
	
	I3Particle muon = primary.Clone();
	muon.SetType(I3Particle::MuMinus);
	muon.SetLocationType(I3Particle::InIce);
	muon.SetShape(I3Particle::Null);
	
	I3MCTreeUtils::AddPrimary(tree, primary);
	I3MCTreeUtils::AppendChild(tree, primary, muon);
}

double
Floodlight::GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	std::pair<double, double> steps = surface_->GetIntersection(axis.GetPos(), axis.GetDir());
	// Bail if the axis doesn't intersect the surface, or there's more than 1 muon.
	if (!std::isfinite(steps.first) || bundle.size() != 1)
		return -std::numeric_limits<double>::infinity();
	
	return energyGenerator_->GetLog(bundle.front().energy) - std::log(surface_->GetAcceptance(-1, 1));
}

}

I3_SERIALIZABLE(I3MuonGun::Floodlight);
