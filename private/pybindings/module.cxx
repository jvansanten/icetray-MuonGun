
#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>

namespace bp = boost::python;

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/physics/I3Particle.h>

void
register_extras()
{
	bp::class_<std::vector<I3Particle::ParticleType> >("I3ParticleTypeSeries")
	    .def(bp::dataclass_suite<std::vector<I3Particle::ParticleType> >())
	;
}

#define REGISTER_THESE_THINGS                                           \
  (histogram)(TrackBinner)(MuonPropagator)(extras)(I3MuonGun)           \
  (CompactTrack)(RadialDistribution)(EnergyDistribution)                \
  (Flux)(Generator)(Surface)(WeightCalculator)(CanCan)

#define I3_REGISTRATION_FN_DECL(r, data, t) void BOOST_PP_CAT(register_,t)();
#define I3_REGISTER(r, data, t) BOOST_PP_CAT(register_,t)();
BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_THESE_THINGS)

I3_PYTHON_MODULE(MuonGun)
{
	load_project("libMuonGun", false);
	BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_THESE_THINGS);
}

