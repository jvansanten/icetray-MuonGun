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

namespace simclasses {

Surface::~Surface() {}

template <typename Archive>
void
Surface::serialize(Archive &ar __attribute__ ((unused)), unsigned version __attribute__ ((unused)))
{}

SamplingSurface::~SamplingSurface() {}

double
SamplingSurface::SampleImpactRay(I3Position &impact, I3Direction &dir, I3RandomService &rng,
    double cosMin, double cosMax) const
{
	dir = SampleDirection(rng, cosMin, cosMax);
	impact = SampleImpactPosition(dir, rng);

	// Calculate projected area
	return GetArea(dir);
}

I3Direction
SamplingSurface::SampleDirection(I3RandomService &rng,
    double cosMin, double cosMax) const
{
	// Sample a direction proportional to the projected area 
	// of the surface.
	double maxarea = GetMaximumArea();
	I3Direction sampled_dir;
	do {
		sampled_dir = I3Direction(acos(rng.Uniform(cosMin, cosMax)),
		    rng.Uniform(0, 2*M_PI));
	} while (rng.Uniform(0, maxarea) > GetArea(sampled_dir));
	
	return sampled_dir;
}

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Surface", base_object<Surface>(*this));
}

// Find the distances to the points of intersection with a centered at (0,0,0)
// and aligned along the z axis. Adapted from:
// http://code.icecube.wisc.edu/svn/projects/mmc/trunk/src/tfa/Amanda.java
// (D. Chirkin)

namespace detail {

template <typename T>
std::pair<double, double>
CylinderBase<T>::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	std::pair<double, double> h(Surface::no_intersection()), r(Surface::no_intersection());
	
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
				h = Surface::no_intersection();
		// Perfectly vertical tracks never intersect the sides
		} else if (sinth == 0) {
			if (hypot(x, y) >= radius_)
				h = Surface::no_intersection();
		// For general tracks, take the last entrace and first exit
		} else {
			if (h.first >= r.second || h.second <= r.first)
				h = Surface::no_intersection();
			else {
				h.first = std::max(r.first, h.first);
				h.second = std::min(r.second, h.second);
			}
		}
	}
	
	return h;
}

template <typename T>
double
CylinderBase<T>::GetArea(const I3Direction &dir) const
{
	return GetAreaForZenith(-dir.GetZ());
}

template <typename T>
double
CylinderBase<T>::GetAreaForZenith(double coszen) const
{
	double cap = M_PI*radius_*radius_;
	double sides = 2*radius_*length_;
	return cap*fabs(coszen) + sides*sqrt(1.-coszen*coszen);
}

template <typename T>
double
CylinderBase<T>::GetMaximumArea() const
{
	double thetaMax = atan(2*length_/(M_PI*radius_));
	return GetAreaForZenith(cos(thetaMax));
}

template <typename T>
I3Direction
CylinderBase<T>::SampleDirection(I3RandomService &rng,
    double cosMin, double cosMax) const
{
	// Sample a direction proportional to the projected area 
	// of the surface.
	double coszen;
	double maxarea = GetMaximumArea();
	do {
		coszen = rng.Uniform(cosMin, cosMax);
	} while (rng.Uniform(0, maxarea) > GetAreaForZenith(coszen));
	
	return I3Direction(acos(coszen), rng.Uniform(0, 2*M_PI));
}

template <typename T>
I3Position
CylinderBase<T>::SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const
{
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
	I3Position impact(y, x, 0);
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
	
	return impact;
}

template <typename Base>
template <typename Archive>
void
CylinderBase<Base>::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Base", base_object<Base>(*this));
	ar & make_nvp("Length", length_);
	ar & make_nvp("Radius", radius_);
	ar & make_nvp("Center", center_);
}

}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

Cylinder::~Cylinder() {}

AxialCylinder::AxialCylinder(double length, double radius, I3Position center)
    : length_(length/2.,length/2.), radius_(radius), center_(center)
{}

AxialCylinder::AxialCylinder(double lengthBefore, double lengthAfter, double radius, I3Position center)
    : length_(lengthBefore,lengthAfter), radius_(radius), center_(center)
{}

std::pair<double, double>
AxialCylinder::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	// Distance to point of closest approach to the center
	double to_center = (center_ - p)*dir;
	// Check distance of closest approach against cylinder radius
	if ((p + to_center*dir - center_).Magnitude() > radius_)
		return no_intersection();
	else
		return std::make_pair(to_center-length_.first, to_center+length_.second);
}

double
AxialCylinder::GetArea(const I3Direction &dir __attribute__((unused))) const
{
	return M_PI*radius_*radius_;
}

double
AxialCylinder::GetAcceptance(double cosMin, double cosMax) const
{
	return M_PI*radius_*radius_*(cosMax-cosMin);
}

double
AxialCylinder::GetMaximumArea() const
{
	return M_PI*radius_*radius_;
}

I3Direction
AxialCylinder::SampleDirection(I3RandomService &rng, double cosMin, double cosMax) const
{
	return I3Direction(std::acos(rng.Uniform(cosMin, cosMax)), rng.Uniform(0, 2*M_PI));
}

I3Position
AxialCylinder::SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const
{
	// Choose a position in a circle in axis-centered coordinates
	I3Position impact(std::sqrt(rng.Uniform(0, radius_*radius_)), 0, 0);
	impact.RotateZ(rng.Uniform(0, 2*M_PI));
	
	// Rotate into the transverse plane
	impact.RotateY(dir.GetZenith());
	impact.RotateZ(dir.GetAzimuth());
	// Shift from cylinder-centered to real coordinates
	impact += center_;
	// Shift back to the entry point
	impact -= length_.first*dir;
	
	return impact;
}

template <typename Archive>
void
AxialCylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("SamplingSurface", base_object<SamplingSurface>(*this));
	ar & make_nvp("Length", length_);
	ar & make_nvp("Radius", radius_);
	ar & make_nvp("Center", center_);
}

AxialCylinder::~AxialCylinder() {}

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

template <typename Archive>
void
Sphere::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Surface", base_object<Surface>(*this));
	ar & make_nvp("OriginDepth", originDepth_);
	ar & make_nvp("Radius", radius_);
}

Sphere::~Sphere() {}

}

namespace I3MuonGun {

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<simclasses::SamplingSurface>(*this));
}

SamplingSurface::~SamplingSurface() {}

std::pair<double, double>
Cylinder::GetZRange() const
{
	return std::make_pair(GetCenter().GetZ() - GetLength()/2., GetCenter().GetZ() + GetLength()/2.);
}

double
Cylinder::GetTopArea() const
{
	return M_PI*GetRadius()*GetRadius();
}

double
Cylinder::GetSideArea() const
{
	return 2*GetRadius()*GetLength();
}

bool
Cylinder::operator==(const SamplingSurface &s) const
{
	const Cylinder *other = dynamic_cast<const Cylinder*>(&s);
	if (!other)
		return false;
	else 
		return (GetRadius() == other->GetRadius() &&
		    GetLength() == other->GetLength() &&
		    GetCenter() == other->GetCenter());
}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 1)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	if (version == 0) {
		// Class hierarchy changed from v0 to v1, so we have to
		// deserialize by hand
		double radius, length;
		I3Position center;
		ar & make_nvp("SamplingSurface", base_object<simclasses::SamplingSurface>(*this));
		ar & make_nvp("Length", length);
		ar & make_nvp("Radius", radius);
		ar & make_nvp("Center", center);
		SetLength(length);
		SetRadius(radius);
		SetCenter(center);
	} else {
		ar & make_nvp("Base", base_object<Base>(*this));
	}
}

}

// explicitly instantiate the base classes used
template class simclasses::detail::CylinderBase<simclasses::SamplingSurface>;
template class simclasses::detail::CylinderBase<I3MuonGun::SamplingSurface>;
template class I3MuonGun::detail::UprightSurface<simclasses::detail::CylinderBase<I3MuonGun::SamplingSurface> >;

I3_SERIALIZABLE(simclasses::Surface);
I3_SERIALIZABLE(simclasses::SamplingSurface);
I3_SERIALIZABLE(simclasses::Cylinder);
I3_SERIALIZABLE(simclasses::Sphere);
I3_SERIALIZABLE(simclasses::AxialCylinder);

I3_SERIALIZABLE(I3MuonGun::SamplingSurface);
I3_SERIALIZABLE(I3MuonGun::Cylinder);


