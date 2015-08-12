/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <simclasses/Cylinder.h>
#include <simclasses/ExtrudedPolygon.h>
#include <simclasses/AxialCylinder.h>
#include <simclasses/Sphere.h>

#include <MuonGun/Cylinder.h>
#include <MuonGun/ExtrudedPolygon.h>

#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <phys-services/I3RandomService.h>
#include <icetray/python/dataclass_suite.hpp>
#include <MuonGun/Flux.h>
#include <icetray/python/gil_holder.hpp>
#include "utils.h"

static double IntegrateFlux(const I3MuonGun::SamplingSurface &s, I3MuonGun::FluxPtr flux, unsigned m, double cosMin, double cosMax)
{
	return s.IntegrateFlux(boost::bind(boost::cref(*flux), _1, _2, m), cosMin, cosMax);
}

static boost::python::tuple SampleImpactRay(const simclasses::SamplingSurface &s, I3RandomServicePtr rng, double cosMin, double cosMax)
{
	I3Position pos;
	I3Direction dir;
	s.SampleImpactRay(pos, dir, *rng, cosMin, cosMax);
	return boost::python::make_tuple(pos, dir);
}

namespace simclasses {

using namespace boost::python;

class PySurface : public Surface, public wrapper<Surface> {
public:
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const
	{
		boost::python::detail::gil_holder lock;
		return get_override("GetIntersection")(p, dir);
	}
	
	virtual bool operator==(const Surface &s) const
	{
		log_fatal("Python Surfaces can't be compared.");
		return false;
	}
};

I3_POINTER_TYPEDEFS(PySurface);

}

void register_base_Surface()
{
	using namespace simclasses;
	using namespace boost::python;
	
	class_<PySurface, PySurfacePtr, boost::noncopyable>("Surface")
	    .def("intersection", &Surface::GetIntersection)
	;
	
	implicitly_convertible<SurfacePtr, SurfaceConstPtr>();
	
	class_<SamplingSurface, SamplingSurfacePtr, bases<Surface>, boost::noncopyable>("SamplingSurface", no_init)
	    .def("area", &SamplingSurface::GetArea, (arg("dir")))
	    .def("maximum_area", &SamplingSurface::GetMaximumArea)
	    .def("sample_impact_ray", &SampleImpactRay, (arg("self"), arg("rng"), arg("cosMin")=0, arg("cosMax")=1))
	    .def("sample_impact_position", &SamplingSurface::SampleImpactPosition, (arg("self"), arg("dir"), arg("rng")))
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
	
	class_<ExtrudedPolygon, ExtrudedPolygonPtr, bases<SamplingSurface> >("ExtrudedPolygon",
	    init<const std::vector<I3Position> &, double>((arg("points"), arg("padding")=0)))
	    .add_property("x", &ExtrudedPolygon::GetX)
	    .add_property("y", &ExtrudedPolygon::GetY)
	    .add_property("z", &ExtrudedPolygon::GetZ)
	;
	
	class_<Sphere, bases<Surface> >("Sphere", init<double, double>())
	;
	
	class_<AxialCylinder, bases<SamplingSurface> >("AxialCylinder",
	    init<double,double,I3Position>((bp::arg("length"), bp::arg("radius"), bp::arg("center")=I3Position(0,0,0))))
	    .def(init<double,double,double,I3Position>((bp::arg("lengthBefore"), "lengthAfter", "radius", bp::arg("center")=I3Position(0,0,0))))
	;
}

void register_Surface()
{
	register_base_Surface();
	
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<SamplingSurface, SamplingSurfacePtr, bases<simclasses::SamplingSurface>, boost::noncopyable>("SamplingSurface", no_init)
	    .def("integrate_flux", &IntegrateFlux, (arg("self"), arg("flux"), arg("m")=1u, arg("cosMin")=0, arg("cosMax")=1))
	    DEF("acceptance", &SamplingSurface::GetAcceptance, (arg("cosMin")=0., arg("cosMax")=1.))
	;
	
	implicitly_convertible<SamplingSurfacePtr, SamplingSurfaceConstPtr>();
	
	class_<Cylinder, CylinderPtr, bases<SamplingSurface> >("Cylinder",
	    init<double, double, I3Position>((arg("length"), arg("radius"), arg("center")=I3Position(0,0,0))))
	    .def(copy_suite<Cylinder>())
	    #define PROPS (Length)(Radius)(Center)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, Cylinder, PROPS)
	    #undef PROPS
	    .def(self == self)
	;

	implicitly_convertible<CylinderPtr, CylinderConstPtr>();

	class_<ExtrudedPolygon, ExtrudedPolygonPtr, bases<SamplingSurface> >("ExtrudedPolygon",
	    init<const std::vector<I3Position> &, double>((arg("points"), arg("padding")=0)))
	    .add_property("x", &ExtrudedPolygon::GetX)
	    .add_property("y", &ExtrudedPolygon::GetY)
	    .add_property("z", &ExtrudedPolygon::GetZ)
	;
}
