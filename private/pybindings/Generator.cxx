/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Generator.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/python/dataclass_suite.hpp>
#include <boost/serialization/list.hpp>

void register_Generator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<GenerationProbability, GenerationProbabilityPtr, boost::noncopyable>("GenerationProbability", no_init)
	    .def("generated_events", &GenerationProbability::GetGeneratedEvents)
	    .add_property("total_events", &GenerationProbability::GetTotalEvents, &GenerationProbability::SetTotalEvents)
	    .def("__add__",  (GenerationProbabilityPtr (*)(GenerationProbabilityPtr, GenerationProbabilityPtr))(&operator+))
	    .def("__mul__",  (GenerationProbabilityPtr (*)(GenerationProbabilityPtr, size_t))(&operator*))
	    .def("__rmul__", (GenerationProbabilityPtr (*)(GenerationProbabilityPtr, size_t))(&operator*))
	    .def("__imul__", (GenerationProbabilityPtr (*)(GenerationProbabilityPtr, size_t))(&operator*=))
	;
	
	class_<Generator, bases<GenerationProbability>, boost::noncopyable>("Generator", no_init)
	;
	
	class_<BundleConfiguration::value_type>("BundleEntry", init<double, double>())
	;
	
	class_<BundleConfiguration, boost::shared_ptr<BundleConfiguration> >("BundleConfiguration")
	    .def(dataclass_suite<BundleConfiguration>())
	;
}
