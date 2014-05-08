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
	return M_PI*radius_*(radius_*fabs(coszen) + (2*length_/M_PI)*sqrt(1-coszen*coszen));
}

static double integrate_area(double a, double b, double radius, double length)
{
	return (M_PI/2)*radius*(radius*(b*b - a*a)
	    + (2*length/M_PI)*(acos(a) - acos(b) -
	      (sqrt(1-a*a)*a) + sqrt(1-b*b)*b));
}

double
Cylinder::GetTotalArea(double cosMin, double cosMax) const
{
	if (cosMin >= 0 && cosMax >= 0)
		return integrate_area(cosMin, cosMax, radius_, length_);
	else if (cosMin < 0 && cosMax <= 0)
		return integrate_area(-cosMax, -cosMin, radius_, length_);
	else if (cosMin < 0 && cosMax > 0)
		return integrate_area(0, -cosMin, radius_, length_)
		    + integrate_area(0, cosMax, radius_, length_);
	else
		log_fatal("Can't deal with zenith range [%.1e, %.1e]", cosMin, cosMax);
	return NAN;
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
	return 2*radius_*sqrt(1-coszen*coszen);
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
		total += 2*M_PI*Integrate(dN_dOmega, cosMin, cosMax, 1e-3, 1e-3);
	}
	
	// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
	{
		f2 dN = boost::bind(flux, boost::bind(GetDepth, _1), _2);
		f2 dOmega = boost::bind(&Cylinder::GetDifferentialSideArea, this, _2);
		f2 dN_dOmega = detail::multiply<2>(dN, dOmega);
		boost::array<double, 2> low = {{center_.GetZ()  - length_/2., cosMin}};
		boost::array<double, 2> high = {{center_.GetZ() + length_/2., cosMax}};
		total += 2*M_PI*Integrate(dN_dOmega, low, high, 1e-3, 1e-3, 10000u);
	}
	
	return total;
}

double
Cylinder::SampleImpactRay(I3Position &impact, I3Direction &dir, I3RandomService &rng,
    double cosMin, double cosMax) const
{
	// Sample a direction proportional to the projected area 
	// of the surface.
	double coszen;
	double maxarea = GetMaxDifferentialArea();
	do {
		coszen = rng.Uniform(cosMin, cosMax);
	} while (rng.Uniform(0, maxarea) > GetDifferentialArea(coszen));
	dir = I3Direction(acos(coszen), rng.Uniform(0, 2*M_PI));
	
	// The projection of a cylinder onto a plane whose
	// normal is inclined by `zenith` w.r.t to the cylinder
	// axis is a rectangle of width 2*r and height h*sin(theta)
	// capped with two half-ellipses of major axis r and
	// minor axis r*cos(theta). Pick a point from a uniform
	// distribution over this area.
	double a = sin(dir.GetZenith())*length_/2.;
	double b = fabs(cos(dir.GetZenith()))*radius_;
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

bool
Cylinder::operator==(const Surface &s) const
{
	const Cylinder *other = dynamic_cast<const Cylinder*>(&s);
	if (!other)
		return false;
	else 
		return (radius_ == other->radius_ &&
		    length_ == other->length_ &&
		    center_ == other->center_);
}

bool
Sphere::operator==(const Surface &s) const
{
	const Sphere *other = dynamic_cast<const Sphere*>(&s);
	if (!other)
		return false;
	else 
		return (radius_ == other->radius_ &&
		    originDepth_ == other->originDepth_);
}

template <typename Archive>
void
Surface::serialize(Archive &ar, unsigned version)
{}

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("Surface", base_object<Surface>(*this));
}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("SamplingSurface", base_object<SamplingSurface>(*this));
	ar & make_nvp("Length", length_);
	ar & make_nvp("Radius", radius_);
	ar & make_nvp("Center", center_);
}

template <typename Archive>
void
Sphere::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("Surface", base_object<Surface>(*this));
	ar & make_nvp("OriginDepth", originDepth_);
	ar & make_nvp("Radius", radius_);
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

I3_SERIALIZABLE(I3MuonGun::Surface);
I3_SERIALIZABLE(I3MuonGun::SamplingSurface);
I3_SERIALIZABLE(I3MuonGun::Cylinder);
I3_SERIALIZABLE(I3MuonGun::Sphere);
