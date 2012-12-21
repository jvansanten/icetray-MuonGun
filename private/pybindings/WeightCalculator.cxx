
#include <MuonGun/WeightCalculator.h>

void register_WeightCalculator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<WeightCalculator>("WeightCalculator", init<SamplingSurfacePtr, FluxPtr,
	    RadialDistributionPtr, EnergyDistributionPtr, GenerationProbabilityPtr>())
	;
}
