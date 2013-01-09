/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/CanCan.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>

void
register_CanCan()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<StaticSurfaceInjector, bases<Generator> >("StaticSurfaceInjector")
		#define PROPS (Surface)(Flux)(RadialDistribution)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, StaticSurfaceInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &StaticSurfaceInjector::GetTotalRate)
	;

#if 1
	class_<EnergyDependentSurfaceInjector, bases<Generator> >("EnergyDependentSurfaceInjector")
		.def("total_rate", &EnergyDependentSurfaceInjector::GetTotalRate)
		.def("get_surface", &EnergyDependentSurfaceInjector::GetSurface)
	;
#endif
}