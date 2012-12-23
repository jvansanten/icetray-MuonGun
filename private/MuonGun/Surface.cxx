/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Surface.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <phys-services/I3RandomService.h>
#include <boost/bind.hpp>

namespace I3MuonGun {

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

inline std::pair<double, double>
no_intersection()
{
	return std::make_pair(NAN, NAN);
}

static const double SurfaceRadius = 6371300+2834;

}

Surface::~Surface() {}

SamplingSurface::~SamplingSurface() {}

// Find the distances to the points of intersection with a centered at (0,0,0)
// and aligned along the z axis. Adapted from:
// http://code.icecube.wisc.edu/svn/projects/mmc/trunk/src/tfa/Amanda.java
// (D. Chirkin)

std::pair<double, double>
Cylinder::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	std::pair<double, double> h(no_intersection()), r(no_intersection());
	
	double x = p.GetX()-center_.GetX();
	double y = p.GetY()-center_.GetY();
	double z = p.GetZ()-center_.GetZ();
	
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
				h = no_intersection();
		// Perfectly vertical tracks never intersect the sides
		} else if (sinth == 0) {
			if (hypot(x, y) >= radius_)
				h = no_intersection();
		// For general tracks, take the last entrace and first exit
		} else {
			if (h.first >= r.second || h.second <= r.first)
				h = no_intersection();
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
	return GetDepth(center_.GetZ() + length_/2.);
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
		f1 dN = boost::bind<double>(flux, GetDepth(center_.GetZ() + length_/2.), _1);
		f1 dOmega = boost::bind(&Cylinder::GetDifferentialTopArea, this, _1);
		f1 dN_dOmega = detail::multiply<1>(dN, dOmega);
		total += 2*M_PI*Integrate(dN_dOmega, cosMin, cosMax, 1e-6, 1e-6);
	}
	
	// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
	{
		f2 dN = boost::bind(flux, boost::bind(GetDepth, _1), _2);
		f2 dOmega = boost::bind(&Cylinder::GetDifferentialSideArea, this, _2);
		f2 dN_dOmega = detail::multiply<2>(dN, dOmega);
		boost::array<double, 2> low = {{center_.GetZ()  - length_/2., 0.}};
		boost::array<double, 2> high = {{center_.GetZ() + length_/2., 1.}};
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
	// Rotate into the transverse plane
	impact = I3Position(y, x, 0);
	impact.RotateY(dir.GetZenith());
	impact.RotateZ(dir.GetAzimuth());
	// Shift from cylinder-centered to real coordinates
	impact.SetX(impact.GetX() + center_.GetX());
	impact.SetY(impact.GetY() + center_.GetY());
	impact.SetZ(impact.GetZ() + center_.GetZ());
	// Project back to the entry point
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
	std::pair<double, double> h(no_intersection());
	
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

}