/** $Id: Surface.h 135924 2015-08-06 14:17:50Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 135924 $
 * $Date: 2015-08-06 16:17:50 +0200 (Thu, 06 Aug 2015) $
 */

#ifndef I3MUONGUN_SURFACE_H_INCLUDED
#define I3MUONGUN_SURFACE_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>
#include <icetray/I3Units.h>
#include <boost/function.hpp>
#include <dataclasses/I3Position.h>
#include <MuonGun/UprightSurface.h>

class I3Direction;
class I3RandomService;

namespace simclasses {

/**
 * @brief A closed surface
 *
 * Surface knows how to find the intersections of a ray with itself.
 */
class Surface {
public:
	virtual ~Surface();
	/**
	 * Find the points where a ray intersects the surface
	 *
	 * @param[in] p   The origin of the ray
	 * @param[in] dir The direction of the ray
	 * @returns a pair of distances from the origin of the ray to the
	 *          intersections with the surface. A distance of NAN means
	 *          that the ray never intersects the surface, and negative
	 *          distances mean that the intersection in question is
	 *          "behind" the origin.
	 */
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const = 0;
	std::pair<double, double> no_intersection() const
	{
		return std::make_pair(NAN, NAN);
	}

private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Surface);

/**
 * @brief A surface that can sample from the distribution of impact positions
 *        and directions of rays from an isotropic flux.
 */
class SamplingSurface : public Surface {
public:
	virtual ~SamplingSurface();
	/**
	 * @brief Calculate the projected area in m^2
	 *
	 * @param[in] dir The direction onto which to project the area
	 */
	virtual double GetArea(const I3Direction &dir) const = 0;
	virtual double GetMaximumArea() const = 0;

	virtual I3Direction SampleDirection(I3RandomService &rng, double cosMin=0, double cosMax=1) const;
	virtual I3Position SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const = 0;

	 /**
	  * Sample an impact point and direction from an isotropic flux
	  *
	  * @param[out] pos    Impact point
	  * @param[out] dir    direction
	  * @param[in]  rng    Random number generator
	  * @param[in]  cosMin cosine of the maximum zenith angle to consider
	  * @param[in]  cosMax cosine of the minimum zenith angle to consider
	  * @returns the projected area along the chosen zenith angle
	  */
	virtual double SampleImpactRay(I3Position &pos, I3Direction &dir, I3RandomService &rng, double cosMin=0, double cosMax=1) const;
private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(SamplingSurface);

namespace detail {

template <typename T>
class CylinderBase : public T {
public:
	CylinderBase(double length, double radius, I3Position center=I3Position(0,0,0)) : length_(length), radius_(radius), center_(center) {};

	// simclasses::Surface interface
	std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;

	// simclasses::SamplingSurface interfact
	double GetArea(const I3Direction &dir) const;
	double GetMaximumArea() const;
	
	// SamplingSurface interface
	I3Direction SampleDirection(I3RandomService &rng, double cosMin=0, double cosMax=1) const;
	I3Position SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const;
	void SetLength(double v) { length_ = v; }
	double GetLength() const { return length_; }
	
	void SetRadius(double v) { radius_ = v; }
	double GetRadius() const { return radius_; }
	
	void SetCenter(const I3Position &v) { center_ = v; }
	I3Position GetCenter() const { return center_; }

protected:
	CylinderBase() {}

private:
	double GetAreaForZenith(double cosZen) const;

	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);

	double length_, radius_;
	I3Position center_;
};

}

/**
 * @brief A cylinder aligned with the z axis
 */
class Cylinder : public detail::CylinderBase<SamplingSurface> {
private:
	typedef detail::CylinderBase<SamplingSurface> Base;
public:
	virtual ~Cylinder();
	Cylinder(double length, double radius, I3Position center=I3Position(0,0,0)) : CylinderBase(length, radius, center) {};

private:
	Cylinder() {}
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Cylinder);

/**
 * @brief A sphere with its origin at a given vertical depth
 */
class Sphere : public Surface {
public:
	virtual ~Sphere();
	Sphere(double originDepth, double radius) : originDepth_(originDepth), radius_(radius) {};
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	
private:
	Sphere() {}
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	double originDepth_, radius_;
};

/**
 * @brief A cylinder aligned with the incoming particle axis a la NuGen
 */
class AxialCylinder : public SamplingSurface {
public:
	virtual ~AxialCylinder();
	AxialCylinder(double length, double radius, I3Position center=I3Position(0,0,0) );
	AxialCylinder(double lengthBefore, double lengthAfter, double radius, I3Position center=I3Position(0,0,0) );
	
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	
	// SamplingSurface interface
	virtual double GetArea(const I3Direction &) const;
	virtual double GetMaximumArea() const;
	virtual double GetAcceptance(double cosMin=0, double cosMax=1) const;
	
	virtual I3Direction SampleDirection(I3RandomService &rng, double cosMin=0, double cosMax=1) const;
	virtual I3Position SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const;
private:
	AxialCylinder() {}
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
	
	std::pair<double, double> length_;
	double radius_;
	I3Position center_;
};

I3_POINTER_TYPEDEFS(AxialCylinder);

}

namespace I3MuonGun {

	/**
	 * @brief A surface upon which muon bundles can be generated
	 *
	 * SamplingSurface knows how to calculate its projected area and integrate
	 * a flux over its surface. It is assumed to be azimuthally symmetric, but
	 * its projected area may vary with zenith angle.
	 */
	class SamplingSurface : public simclasses::SamplingSurface {
	public:
		virtual ~SamplingSurface();
		/** Get the integral of area*solid_angle over the given cos(theta) range */
		virtual double GetAcceptance(double cosMin=0, double cosMax=1) const = 0;
		/** Get the minimum vertical depth the surface occupies */
		virtual double GetMinDepth() const = 0;
	
		/** 
		 * Integrate a flux (defined in terms of dN/dOmega(depth [km], cos(theta)))
		 * over the outer surface
		 */
		virtual double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const = 0;
	
		virtual bool operator==(const SamplingSurface&) const = 0;
	private:
		friend class boost::serialization::access;
		template <typename Archive>
		void serialize(Archive &, unsigned);
	};
	
	I3_POINTER_TYPEDEFS(SamplingSurface);

	/**
	 * @brief A cylinder aligned with the z axis
	 */
	class Cylinder : public detail::UprightSurface<simclasses::detail::CylinderBase<SamplingSurface > > {
	private:
		typedef simclasses::detail::CylinderBase<SamplingSurface > CylinderBase;
		typedef detail::UprightSurface<CylinderBase > Base;
	public:
		Cylinder(double length, double radius, I3Position center=I3Position(0,0,0)) : Base(length, radius, center)
		{}
	
		// SamplingSurface interface
		bool operator==(const SamplingSurface&) const;
	
		double GetLength() const { return CylinderBase::GetLength(); };
	
	protected:
		// UprightSurface interface
		double GetTopArea() const;
		double GetSideArea() const;
		std::pair<double, double> GetZRange() const;
		
	private:
		Cylinder() {}
		
		friend class boost::serialization::access;
		template <typename Archive>
		void serialize(Archive &, unsigned);
	};

	I3_POINTER_TYPEDEFS(Cylinder);
}

BOOST_CLASS_VERSION(simclasses::Surface, 0);
BOOST_CLASS_VERSION(simclasses::SamplingSurface, 0);
BOOST_CLASS_VERSION(simclasses::Cylinder, 0);
BOOST_CLASS_VERSION(simclasses::Sphere, 0);
BOOST_CLASS_VERSION(simclasses::AxialCylinder, 0);

BOOST_CLASS_VERSION(I3MuonGun::SamplingSurface, 0);
BOOST_CLASS_VERSION(I3MuonGun::Cylinder, 1);

#endif
