
#include <MuonGun/RadialDistribution.h>

void
register_RadialDistribution()
{
	namespace bp = boost::python;
	using namespace I3MuonGun;
	
	bp::class_<Distribution, DistributionPtr,
	    bp::bases<I3FrameObject> >("Distribution")
	    .add_property("random_service", &Distribution::GetRandomService, &Distribution::SetRandomService)
	;
	
	bp::class_<Sample>("Sample", bp::init<double,double>())
	    .def_readwrite("value", &Sample::value)
	    .def_readwrite("prob", &Sample::prob)
	;
	
	bp::class_<RadialDistribution, RadialDistributionPtr,
	    bp::bases<Distribution>, boost::noncopyable>("RadialDistribution", bp::no_init)
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
