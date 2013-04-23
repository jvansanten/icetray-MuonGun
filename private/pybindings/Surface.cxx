/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Surface.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <icetray/python/dataclass_suite.hpp>
#include <MuonGun/Flux.h>

static double IntegrateFlux(const I3MuonGun::SamplingSurface &s, I3MuonGun::FluxPtr flux, unsigned m, double cosMin, double cosMax)
{
	return s.IntegrateFlux(boost::bind(boost::cref(*flux), _1, _2, m), cosMin, cosMax);
}

void register_Surface()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<Surface, SurfacePtr, boost::noncopyable>("Surface", no_init)
	    .def("intersection", &Surface::GetIntersection)
	;
	
	implicitly_convertible<SurfacePtr, SurfaceConstPtr>();
	
	class_<SamplingSurface, SamplingSurfacePtr, bases<Surface>, boost::noncopyable>("SamplingSurface", no_init)
	    .def("differential_area", &SamplingSurface::GetDifferentialArea)
	    .def("total_area", &SamplingSurface::GetTotalArea, (arg("cosMin")=0., arg("cosMax")=1.))
	    .def("integrate_flux", &IntegrateFlux, (arg("self"), arg("flux"), arg("m")=1u, arg("cosMin")=0, arg("cosMax")=1))
	;
	
	implicitly_convertible<SamplingSurfacePtr, SamplingSurfaceConstPtr>();
	
	class_<Cylinder, CylinderPtr, bases<SamplingSurface> >("Cylinder",
	    init<double, double, I3Position>((arg("length"), arg("radius"), arg("center")=I3Position(0,0,0))))
	    .def(copy_suite<Cylinder>())
	    #define PROPS (Length)(Radius)(Center)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, Cylinder, PROPS)
	    #undef PROPS
	;
	
	implicitly_convertible<CylinderPtr, CylinderConstPtr>();
	
	class_<Sphere, bases<Surface> >("Sphere", init<double, double>())
	;
}
