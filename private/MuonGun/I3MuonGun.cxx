
#include <MuonGun/I3MuonGun.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3Logging.h>

#include <gsl/gsl_integration.h>
#include <cubature/cubature.h>

namespace I3MuonGun {

// SET_LOGGER("I3MuonGun");

namespace {

inline void
sort(std::pair<double, double> &pair)
{
	if (pair.first > pair.second) {
		double aux = pair.first;
		pair.first = pair.second;
		pair.second = aux;
	}
}

}

// Find the distances to the points of intersection with a centered at (0,0,0)
// and aligned along the z axis. Adapted from:
// http://code.icecube.wisc.edu/svn/projects/mmc/trunk/src/tfa/Amanda.java
// (D. Chirkin)
std::pair<double, double>
CylinderIntersection(const I3Position &p, const I3Direction &dir, double length, double radius)
{
	std::pair<double, double> h(0, 0), r(0, 0);
	
	double x = p.GetX();
	double y = p.GetY();
	double z = p.GetZ();
	
	double sinph = sin(dir.GetAzimuth());
	double cosph = cos(dir.GetAzimuth());
	double sinth = sin(dir.GetZenith());
	double costh = cos(dir.GetZenith());
	
	double b = x*cosph + y*sinph;
	double d = b*b + radius*radius - x*x - y*y;
	
	if (d > 0) {
		d = sqrt(d);
		// down-track distance to the endcaps
		if (costh != 0) {
			h.first  = (z - length/2)/costh;
			h.second = (z + length/2)/costh;
			sort(h);
		}
		// down-track distance to the side surfaces
		if (sinth != 0) {
			r.first  = (b - d)/sinth;
			r.second = (b + d)/sinth;
			sort(r);
		}
		// Perfectly horizontal tracks never intersect the endcaps
		if (costh == 0) {
			if ((z > -length/2) && (z < length/2))
				h = r;
			else
				h = std::make_pair(0, 0);
		// Perfectly vertical tracks never intersect the sides
		} else if (sinth == 0) {
			if (hypot(x, y) >= radius)
				h = std::make_pair(0, 0);
		// For general tracks, take the last entrace and first exit
		} else {
			if (h.first >= r.second || h.second <= r.first)
				h = std::make_pair(0, 0);
			else {
				h.first = std::max(r.first, h.first);
				h.second = std::min(r.second, h.second);
			}
		}
	}
	
	return h;
}

namespace {

double gsl_thunk(double x, void *p)
{
	typedef boost::function<double (double)> func_t;
	func_t *f = static_cast<func_t*>(p);
	return (*f)(x);
}

void integrate_thunk(unsigned ndims, const double *x, void *p, unsigned fdims, double *fval)
{
	assert(ndims == 2);
	assert(fdims == 1);
	
	typedef boost::function<double (double, double)> func_t;
	func_t *f = static_cast<func_t*>(p);
	fval[0] = (*f)(x[0], x[1]);
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

double Integrate(boost::function<double (double, double)> f, double low[2], double high[2], double epsabs, double epsrel, size_t maxcall)
{
	double result, error;
	// Special case: 2-dimensional, scalar-valued.
	unsigned fdims = 1;
	unsigned ndims = 2;
	
	int err = adapt_integrate(fdims, &integrate_thunk, &f, ndims, low, high, maxcall, epsabs, epsrel, &result, &error);
	
	return result;
}


}
