/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/CanCan.h>
#include <MuonGun/Floodlight.h>
#include <MuonGun/EnergyDependentSurfaceInjector.h>
#include <icetray/python/function.hpp>

void
register_CanCan()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<StaticSurfaceInjector, bases<Generator> >("StaticSurfaceInjector")
		.def(init<SamplingSurfacePtr, FluxPtr,
		    boost::shared_ptr<OffsetPowerLaw>, RadialDistributionPtr>())
		#define PROPS (Flux)(RadialDistribution)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, StaticSurfaceInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &StaticSurfaceInjector::GetTotalRate)
	;

#if 1
	def_function<boost::function<SamplingSurfacePtr (double)> >("SurfaceScalingFunction");
	class_<EnergyDependentSurfaceInjector, bases<Generator> >("EnergyDependentSurfaceInjector",
	    init<FluxPtr, RadialDistributionPtr, boost::shared_ptr<OffsetPowerLaw>, boost::function<SamplingSurfacePtr (double)> >
	    ((arg("flux")=FluxPtr(), arg("radius")=RadialDistributionPtr(),
	    arg("energy")=boost::shared_ptr<OffsetPowerLaw>(), arg("scaling")=boost::function<SamplingSurfacePtr (double)>())))
		.def("total_rate", &EnergyDependentSurfaceInjector::GetTotalRate)
		.def("target_surface", &EnergyDependentSurfaceInjector::GetTargetSurface)
		#define PROPS (Scaling)(Flux)(EnergyDistribution)(RadialDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, EnergyDependentSurfaceInjector, PROPS)
		#undef PROPS
	;
#endif
	class_<Floodlight, boost::shared_ptr<Floodlight>, bases<Generator> >("Floodlight")
	;
}
