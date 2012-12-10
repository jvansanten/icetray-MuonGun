
#include <MuonGun/CanCan.h>

void
register_CanCan()
{
	namespace bp = boost::python;
	using namespace I3MuonGun;
	
	// bp::class_<CanCan, boost::shared_ptr<CanCan>,
	//     bp::bases<Distribution>, boost::noncopyable>("CanCan", bp::init<double, double>(
	//     (bp::arg("radius"), bp::arg("length"))))
	//     .def("generate_bundle", &CanCan::GenerateBundle)
	//     .def("get_depth", &CanCan::GetDepth)
	//     .add_property("total_rate", &CanCan::GetTotalRate)
	// 	    
	// ;
}