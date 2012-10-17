
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
	
	bp::class_<RadialDistribution, RadialDistributionPtr,
	    bp::bases<Distribution>, boost::noncopyable>("RadialDistribution", bp::no_init)
	    .def("generation_probability", &RadialDistribution::GetGenerationProbability)
	    .def("generate", &RadialDistribution::Generate)
	;
	
	bp::class_<BMSSRadialDistribution, BMSSRadialDistributionPtr,
	    bp::bases<RadialDistribution> >("BMSSRadialDistribution")
	;
}
