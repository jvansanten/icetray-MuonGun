/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/TrackBinner.h>

using namespace I3MuonGun;
using namespace I3MuonGun::histogram;
namespace bp = boost::python;

static boost::shared_ptr<histogram_base> get_primary(TrackBinner &t) { return t.primary_; }
static boost::shared_ptr<histogram_base> get_multiplicity(TrackBinner &t) { return t.multiplicity_; }
static boost::shared_ptr<histogram_base> get_radius(TrackBinner &t) { return t.radius_; }
static boost::shared_ptr<histogram_base> get_energy(TrackBinner &t) { return t.energy_; }

static boost::shared_ptr<histogram_base> get_nue(NeutrinoBinner &t) { return t.nu_e_; }
static boost::shared_ptr<histogram_base> get_numu(NeutrinoBinner &t) { return t.nu_mu_; }

void register_TrackBinner()
{
	bp::class_<TrackBinner>("TrackBinner", bp::init<double,double,unsigned>((
	    bp::arg("mindepth")=1.0, bp::arg("maxdepth")=5.0, bp::arg("steps")=9)))
	    .def("consume", &TrackBinner::Consume)
	    .add_property("primary", &get_primary)
	    .add_property("multiplicity", &get_multiplicity)
	    .add_property("radius", &get_radius)
	    .add_property("energy", &get_energy)
	;
	
	bp::class_<NeutrinoBinner>("NeutrinoBinner")
	    .def("consume", &NeutrinoBinner::Consume)
	    .add_property("nu_e", &get_nue)
	    .add_property("nu_mu", &get_numu)
	;
}