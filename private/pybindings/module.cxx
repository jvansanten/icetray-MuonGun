
#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>

namespace bp = boost::python;

#define REGISTER_THESE_THINGS                                           \
  (histogram)(TrackBinner)(MuonPropagator)(I3MuonGun)                   \
  (CompactTrack)(RadialDistribution)(EnergyDistribution)                \
  (Flux)(Generator)(Surface)(WeightCalculator)(CanCan)                  \
  (CORSIKAGenerationProbability)

#define I3_REGISTRATION_FN_DECL(r, data, t) void BOOST_PP_CAT(register_,t)();
#define I3_REGISTER(r, data, t) BOOST_PP_CAT(register_,t)();
BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_THESE_THINGS)

I3_PYTHON_MODULE(MuonGun)
{
	load_project("libMuonGun", false);
	BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_THESE_THINGS);
}

