/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/WeightCalculator.h>

void register_WeightCalculator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	def("muons_at_surface", &GetMuonsAtSurface);
	
	class_<WeightCalculator>("WeightCalculator", init<SamplingSurfacePtr, FluxPtr,
	    RadialDistributionPtr, EnergyDistributionPtr, GenerationProbabilityPtr>())
	;
}
