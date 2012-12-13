
#include <MuonGun/Generator.h>

void register_Generator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<Generator, boost::noncopyable>("Generator", no_init)
	    .add_property("total_events", &Generator::GetTotalEvents, &Generator::SetTotalEvents)
	;
}
