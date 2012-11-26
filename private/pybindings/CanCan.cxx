
#include <MuonGun/CanCan.h>

void
register_CanCan()
{
	namespace bp = boost::python;
	using namespace I3MuonGun;
	
	bp::class_<CanCan, boost::shared_ptr<CanCan>,
	    bp::bases<Distribution>, boost::noncopyable>("CanCan", bp::init<double, double>())
	    // .def("generation_probability", &RadialDistribution::GetGenerationProbability)
	    .def("generate", &CanCan::Generate)
	    .def("get_depth", &CanCan::GetDepth)
	    .add_property("livetime", &CanCan::GetLivetime)
		    
	;
}