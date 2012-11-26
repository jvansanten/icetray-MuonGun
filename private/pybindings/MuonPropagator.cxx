
#include "MuonGun/MuonPropagator.h"

void
register_MuonPropagator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<MuonPropagator, boost::shared_ptr<MuonPropagator> >("MuonPropagator",
	    init<const std::string&, double, double, double>((
	    arg("medium")="ice", arg("ecut")=-1, arg("vcut")=-1, arg("rho")=1.0)))
	    .def("propagate", &MuonPropagator::propagate)
	    .def("set_seed", &MuonPropagator::SetSeed)
	;
}
