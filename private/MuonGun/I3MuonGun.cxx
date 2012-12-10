
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

inline void
sort(std::pair<double, double> &pair)
{
	if (pair.first > pair.second) {
		double aux = pair.first;
		pair.first = pair.second;
		pair.second = aux;
	}
}

static const double SurfaceRadius = 6371300+2834;

}

// Find the distances to the points of intersection with a centered at (0,0,0)
// and aligned along the z axis. Adapted from:
// http://code.icecube.wisc.edu/svn/projects/mmc/trunk/src/tfa/Amanda.java
// (D. Chirkin)

Surface::~Surface() {}

std::pair<double, double>
Cylinder::GetIntersection(const I3Position &p, const I3Direction &dir) const
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
	double d = b*b + radius_*radius_ - x*x - y*y;
	
	if (d > 0) {
		d = sqrt(d);
		// down-track distance to the endcaps
		if (costh != 0) {
			h.first  = (z - length_/2)/costh;
			h.second = (z + length_/2)/costh;
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
			if ((z > -length_/2) && (z < length_/2))
				h = r;
			else
				h = std::make_pair(0, 0);
		// Perfectly vertical tracks never intersect the sides
		} else if (sinth == 0) {
			if (hypot(x, y) >= radius_)
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

double
Cylinder::GetDifferentialArea(double coszen) const
{
	return M_PI*radius_*(radius_*coszen + (2*length_/M_PI)*sqrt(1-coszen*coszen));
}

double
Cylinder::GetMaxDifferentialArea() const
{
	double thetaMax = atan(2*length_/(M_PI*radius_));
	return GetDifferentialArea(cos(thetaMax));
}

double
Cylinder::GetMinDepth() const
{
	return GetDepth(length_/2.);
}

// dAd Omega/dcos(theta) dphi (only one depth)
double
Cylinder::GetDifferentialTopArea(double coszen) const
{
	return M_PI*radius_*(radius_*coszen);
}

// dAd Omega/dcos(theta) dphi dz (differential also in depth)
double
Cylinder::GetDifferentialSideArea(double coszen) const
{
	return M_PI*radius_*((2/M_PI)*sqrt(1-coszen*coszen));
}

double
Cylinder::IntegrateFlux(boost::function<double (double, double)> flux,
    double cosMin, double cosMax) const
{
        typedef boost::function<double (double)> f1;
	typedef boost::function<double (double, double)> f2;
	
	double total = 0;
	
	// First, integrate to find dN/dt on the cap(s)
	{
		f1 dN = boost::bind<double>(flux, GetDepth(length_/2.), _1);
		f1 dOmega = boost::bind(&Cylinder::GetDifferentialTopArea, this, _1);
		f1 dN_dOmega = detail::multiply<1>(dN, dOmega);
		total += 2*M_PI*Integrate(dN_dOmega, cosMin, cosMax);
	}
	
	// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
	{
		f2 dN = boost::bind(flux, boost::bind(GetDepth, _1), _2);
		f2 dOmega = boost::bind(&Cylinder::GetDifferentialSideArea, this, _2);
		f2 dN_dOmega = detail::multiply<2>(dN, dOmega);
		boost::array<double, 2> low = {{-length_/2., 0.}};
		boost::array<double, 2> high = {{length_/2., 1.}};
		total += 2*M_PI*Integrate(dN_dOmega, low, high, 2e-8, 2e-8, 10000u);
	}
	
	return total;
}

double
Cylinder::SampleImpactRay(I3Position &impact, I3Direction &dir, I3RandomService &rng,
    double cosMin, double cosMax) const
{
	double coszen = rng.Uniform(cosMin, cosMax);
	dir = I3Direction(acos(coszen), rng.Uniform(0, 2*M_PI));
	
	// The projection of a cylinder onto a plane whose
	// normal is inclined by `zenith` w.r.t to the cylinder
	// axis is a rectangle of width 2*r and height h*sin(theta)
	// capped with two half-ellipses of major axis r and
	// minor axis r*cos(theta). Pick a point from a uniform
	// distribution over this area.
	double a = sin(dir.GetZenith())*length_/2.;
	double b = cos(dir.GetZenith())*radius_;
	double x, y;
	do {
		x = radius_*rng.Uniform(-1, 1);
		y = (a + b)*rng.Uniform(-1, 1);
	} while (fabs(y) > a + b*sqrt(1 - (x*x)/(radius_*radius_)));
	impact = I3Position(y, x, 0);
	impact.RotateY(dir.GetZenith());
	impact.RotateZ(dir.GetAzimuth());
	
	// Now, project back to the entry point
	double l = GetIntersection(impact, dir).first;
	impact.SetX(impact.GetX() + l*dir.GetX());
	impact.SetY(impact.GetY() + l*dir.GetY());
	impact.SetZ(impact.GetZ() + l*dir.GetZ());
	
	// Calculate d(A_projected)/d(cos(theta))
	return 2*M_PI*GetDifferentialArea(coszen);
}


std::pair<double, double>
Sphere::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	std::pair<double, double> h(0, 0);
	
	double x = p.GetX();
	double y = p.GetY();
	double z = p.GetZ() - originDepth_;
	
	double sinph = sin(dir.GetAzimuth());
	double cosph = cos(dir.GetAzimuth());
	double sinth = sin(dir.GetZenith());
	double costh = cos(dir.GetZenith());
	
	double b = (x*cosph + y*sinph)*sinth + (z + radius_)*costh;
	double d = b*b - (x*x + y*y + z*(z + 2*radius_));
	
	if (d > 0) {
		d = sqrt(d);
		h.first = b - d;
		h.second = b + d;
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
