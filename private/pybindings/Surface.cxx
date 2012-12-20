
#include <MuonGun/Surface.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

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
	;
	
	implicitly_convertible<SamplingSurfacePtr, SamplingSurfaceConstPtr>();
	
	class_<Cylinder, CylinderPtr, bases<SamplingSurface> >("Cylinder", init<double, double>())
	;
	
	class_<Sphere, bases<Surface> >("Sphere", init<double, double>())
	;
}
