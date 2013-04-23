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
		.def(init<CylinderPtr, FluxPtr,
		    boost::shared_ptr<OffsetPowerLaw>, RadialDistributionPtr>())
		#define PROPS (Flux)(RadialDistribution)(EnergyDistribution)
		BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, StaticSurfaceInjector, PROPS)
		#undef PROPS
		.add_property("total_rate", &StaticSurfaceInjector::GetTotalRate)
	;

	class_<SurfaceScalingFunction, SurfaceScalingFunctionPtr, boost::noncopyable>("SurfaceScalingFunction", no_init)
		.def("__call__", &SurfaceScalingFunction::GetSurface)
	;
	
	class_<BasicSurfaceScalingFunction, BasicSurfaceScalingFunctionPtr, bases<SurfaceScalingFunction> >("BasicSurfaceScalingFunction")
	;

#if 1
	class_<EnergyDependentSurfaceInjector, bases<StaticSurfaceInjector> >("EnergyDependentSurfaceInjector",
	    init<CylinderPtr, FluxPtr, boost::shared_ptr<OffsetPowerLaw>, RadialDistributionPtr, SurfaceScalingFunctionPtr>
	    ((arg("surface")=CylinderPtr(), arg("flux")=FluxPtr(), arg("energy")=boost::shared_ptr<OffsetPowerLaw>(),
	    arg("radius")=RadialDistributionPtr(), arg("scaling")=boost::make_shared<BasicSurfaceScalingFunction>())))
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
