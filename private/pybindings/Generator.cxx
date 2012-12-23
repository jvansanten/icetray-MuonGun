/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Generator.h>
#include <dataclasses/physics/I3Particle.h>

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
}
