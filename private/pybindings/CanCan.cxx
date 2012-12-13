
#include <MuonGun/CanCan.h>

void
register_CanCan()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<StaticSurfaceInjector, bases<Generator, Distribution> >("StaticSurfaceInjector")
		#define PROPS (Surface)(Flux)(RadialDistribution)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, StaticSurfaceInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &StaticSurfaceInjector::GetTotalRate)
	;
}