/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/WeightCalculator.h>
#include <MuonGun/Surface.h>

#include <tableio/converter/pybindings.h>

void register_WeightCalculator()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	def("muons_at_surface", &GetMuonsAtSurface);
	
	class_<WeightCalculator>("WeightCalculator", init<SamplingSurfacePtr, FluxPtr,
	    RadialDistributionPtr, EnergyDistributionPtr, GenerationProbabilityPtr>())
	    .def("__call__", &WeightCalculator::GetWeight)
	    #define PROPS (Surface)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, WeightCalculator, PROPS)
	    #undef PROPS
	;
	
	I3CONVERTER_NAMESPACE(MuonGun);
	I3CONVERTER_EXPORT(MuonBundleConverter, "foo")
	    .def(init<uint32_t, SamplingSurfacePtr>((
	    arg("maxMultiplicity")=25,
	    arg("surface")=boost::make_shared<Cylinder>(1600, 800))))
	;
}
