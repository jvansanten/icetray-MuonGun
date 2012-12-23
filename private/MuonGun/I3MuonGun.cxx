/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/I3MuonGun.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Constants.h>
#include <phys-services/I3RandomService.h>
#include <icetray/I3Logging.h>
#include <icetray/I3Units.h>

#include <gsl/gsl_integration.h>
#include <cubature/cubature.h>

namespace I3MuonGun {
	
double
GetDepth(double z)
{
	return (I3Constants::SurfaceElev - I3Constants::OriginElev - z)/I3Units::km;
}

// SET_LOGGER("I3MuonGun");

namespace {

double gsl_thunk(double x, void *p)
{
	typedef boost::function<double (double)> func_t;
	func_t *f = static_cast<func_t*>(p);
	return (*f)(x);
}

}

double Integrate(boost::function<double (double)> f, double low, double high, double epsabs, double epsrel, size_t limit)
{
	assert(std::isfinite(low));
	assert(std::isfinite(high));
	
	gsl_function gf;
	gf.function = &gsl_thunk;
	gf.params = &f;
	
	double result;
	double abserr;
	
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(limit);
	// Adaptive Gauss-Kronrod 21-point integration rule
	int err = gsl_integration_qags(&gf, low, high, epsabs, epsrel, limit, workspace, &result, &abserr);
	gsl_integration_workspace_free(workspace);
	
	return result;
}

}
