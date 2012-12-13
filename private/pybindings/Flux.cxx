
#include <MuonGun/Flux.h>

void register_Flux()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<Flux, boost::noncopyable>("Flux", no_init)
	    .def("__call__", &Flux::operator(), (args("depth"), "cos_theta", "multiplicity"))
	;
	
	class_<SplineFlux, bases<Flux> >("SplineFlux", init<const std::string&, const std::string&>())
	;
	
	class_<BMSSFlux, bases<Flux> >("BMSSFlux")
	;
}
