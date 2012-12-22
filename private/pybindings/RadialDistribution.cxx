
#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>

void
register_RadialDistribution()
{
	namespace bp = boost::python;
	using namespace I3MuonGun;
	
	bp::class_<RadialDistribution, RadialDistributionPtr,
	    boost::noncopyable>("RadialDistribution", bp::no_init)
	    .def("__call__", &RadialDistribution::operator())
	    .def("generate", &RadialDistribution::Generate)
	;
	
	bp::class_<BMSSRadialDistribution, BMSSRadialDistributionPtr,
	    bp::bases<RadialDistribution> >("BMSSRadialDistribution")
	;
	
	bp::class_<SplineRadialDistribution,
	    bp::bases<RadialDistribution> >("SplineRadialDistribution", bp::init<const std::string &>())
	;
}
