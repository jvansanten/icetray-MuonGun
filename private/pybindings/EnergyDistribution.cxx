/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/EnergyDistribution.h>
#include <MuonGun/RadialDistribution.h>
#include <phys-services/I3RandomService.h>
#include "utils.h"

#include <icetray/python/gil_holder.hpp>

namespace I3MuonGun {

using namespace boost::python;

class PyEnergyDistribution : public EnergyDistribution, public wrapper<EnergyDistribution> {
public:
	virtual double GetLog(double depth, double cos_theta,
	    unsigned multiplicity, double radius, double energy) const
	{
		detail::gil_holder lock;
		return get_override("GetLog")(depth, cos_theta, multiplicity, radius, energy);
	}
	virtual std::pair<double,double> Generate(I3RandomService &rng,
	    double depth, double cos_theta,
	    unsigned multiplicity) const
	{
		detail::gil_holder lock;
		return get_override("Generate")(rng, depth, cos_theta, multiplicity);
	}
	virtual bool operator==(const EnergyDistribution&) const
	{
		return false;
	}
	

};

}

void register_EnergyDistribution()
{
	using namespace boost::python;
	using namespace I3MuonGun;
	
	class_<EnergyDistribution, EnergyDistributionPtr, boost::noncopyable>("EnergyDistribution", no_init)
	    DEF("__call__", &EnergyDistribution::operator(), (arg("depth"), "cos_theta", "multiplicity", "radius", "energy"))
	    .def("generate", (std::pair<double,double> (EnergyDistribution::*)(I3RandomService&,double,double,unsigned) const)&EnergyDistribution::Generate, (arg("rng"), arg("depth"), "cos_theta", "multiplicity"))
	#define PROPS (Min)(Max)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, EnergyDistribution, PROPS)
	#undef PROPS
		// .def("generate", (std::pair<double,double> (EnergyDistribution::*)(I3RandomService&,const RadialDistribution&,double,double,unsigned) const)&EnergyDistribution::Generate, (arg("rng"), arg("radial_dist"), arg("depth"), "cos_theta", "multiplicity"))
	;
	
	class_<SplineEnergyDistribution, boost::shared_ptr<SplineEnergyDistribution>,
	    bases<EnergyDistribution> >("SplineEnergyDistribution",
	    init<const std::string&, const std::string&>((arg("singles"), "bundles")))
	;
	
	class_<BMSSEnergyDistribution, boost::shared_ptr<BMSSEnergyDistribution>,
	    bases<EnergyDistribution> >("BMSSEnergyDistribution")
	    .def("get_spectrum", &BMSSEnergyDistribution::GetSpectrum, (arg("depth"), arg("cos_theta"), arg("multiplicity"), arg("radius")))
	;
	
	// class_<PyEnergyDistribution, boost::noncopyable>("EnergyDistributionBase")
	//     	    DEF("__call__", &EnergyDistribution::operator(), (arg("depth"), "cos_theta", "multiplicity", "radius", "energy"))
	//     	    .def("generate", &EnergyDistribution::Generate, (arg("rng"), arg("depth"), "cos_theta", "multiplicity", "radius"))
	// 	    ;
	
	class_<OffsetPowerLaw, boost::shared_ptr<OffsetPowerLaw> >("OffsetPowerLaw",
	    init<double,double,double,double>((arg("gamma"), "offset", "min", "max")))
	    DEF("__call__", &OffsetPowerLaw::operator(), (arg("energy")))
	    .def("generate", &OffsetPowerLaw::Generate)
	;
}
