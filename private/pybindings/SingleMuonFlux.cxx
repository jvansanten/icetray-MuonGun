
#include <MuonGun/SingleMuonFlux.h>

void register_SingleMuonFlux()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<SingleMuonFlux, bases<Distribution>, boost::noncopyable>("SingleMuonFlux", no_init)
	    .def("__call__", &SingleMuonFlux::operator())
	;
	
	class_<BMSSSingleMuonFlux, bases<SingleMuonFlux> >("BMSSSingleMuonFlux")
	;
	
	class_<AdHocSingleMuonFlux, bases<SingleMuonFlux> >("AdHocSingleMuonFlux")
	;
}
