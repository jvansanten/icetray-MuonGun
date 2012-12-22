
#include <MuonGun/EnergyDistribution.h>
#include <phys-services/I3RandomService.h>

void register_EnergyDistribution()
{
	using namespace boost::python;
	using namespace I3MuonGun;
	
	class_<EnergyDistribution, EnergyDistributionPtr, boost::noncopyable>("EnergyDistribution", no_init)
	    .def("__call__", &EnergyDistribution::operator(), (arg("depth"), "cos_theta", "radius", "multiplicity", "energy"))
	    .def("generate", &EnergyDistribution::Generate, (arg("rng"), arg("depth"), "cos_theta", "multiplicity", "radius"))
	;
	
	class_<SplineEnergyDistribution, boost::shared_ptr<SplineEnergyDistribution>,
	    bases<EnergyDistribution> >("SplineEnergyDistribution",
	    init<const std::string&, const std::string&>((arg("singles"), "bundles")))
	;
	
	class_<OffsetPowerLaw, boost::shared_ptr<OffsetPowerLaw> >("OffsetPowerLaw",
	    init<double,double,double,double>((arg("gamma"), "offset", "min", "max")))
	    .def("__call__", &OffsetPowerLaw::operator())
	    .def("generate", &OffsetPowerLaw::Generate)
	;
}
