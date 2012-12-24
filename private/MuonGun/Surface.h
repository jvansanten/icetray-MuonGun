/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_SURFACE_H_INCLUDED
#define I3MUONGUN_SURFACE_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>
#include <boost/function.hpp>
#include <dataclasses/I3Position.h>

class I3Direction;
class I3RandomService;

namespace I3MuonGun {

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
	};
	
	I3_POINTER_TYPEDEFS(Surface);
	
	/**
	 * @brief A surface upon which muon bundles can be generated
	 *
	 * SamplingSurface knows how to calculate its projected area and integrate
	 * a flux over its surface. It is assumed to be azimuthally symmetric, but
	 * its projected area may vary with zenith angle.
	 */
	class SamplingSurface : public Surface {
	public:
		virtual ~SamplingSurface();
		/** Get d(A_projected)/d(cos(theta)) at a particular cos(theta) */
		virtual double GetDifferentialArea(double coszen) const = 0;
		/** Get the maximum of d(A_projected)/d(cos(theta)) */
		virtual double GetMaxDifferentialArea() const = 0;
		/** Get the minimum vertical depth the surface occupies */
		virtual double GetMinDepth() const = 0;
	
		/** 
		 * Integrate a flux (defined in terms of dN/dOmega(depth [km], cos(theta)))
		 * over the outer surface
		 */
		virtual double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const = 0;
	
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
		virtual double SampleImpactRay(I3Position &pos, I3Direction &dir, I3RandomService &rng, double cosMin=0, double cosMax=1) const = 0;
	};
	
	I3_POINTER_TYPEDEFS(SamplingSurface);

	/**
	 * @brief A cylinder aligned with the z axis
	 */
	class Cylinder : public SamplingSurface {
	public:
		Cylinder(double length, double radius, I3Position center=I3Position(0,0,0)) : length_(length), radius_(radius), center_(center) {};
		
		// Surface interface
		std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
		
		// SamplingSurface interface
		double GetDifferentialArea(double coszen) const;
		double GetMaxDifferentialArea() const;
		double GetMinDepth() const;
		double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const;
		double SampleImpactRay(I3Position &pos, I3Direction &dir, I3RandomService &rng, double cosMin=0, double cosMax=1) const;
		
		void SetLength(double v) { length_ = v; }
		double GetLength() const { return length_; }
		
		void SetRadius(double v) { radius_ = v; }
		double GetRadius() const { return radius_; }
		
		void SetCenter(const I3Position &v) { center_ = v; }
		const I3Position& GetCenter() const { return center_; }
		
	private:
		double GetDifferentialTopArea(double cos_zenith) const;
		double GetDifferentialSideArea(double cos_zenith) const;
	
		double length_, radius_;
		I3Position center_;
	};

	I3_POINTER_TYPEDEFS(Cylinder);

	/**
	 * @brief A sphere with its origin at a given vertical depth
	 */
	class Sphere : public Surface {
	public:
		Sphere(double originDepth, double radius) : originDepth_(originDepth), radius_(radius) {};
		std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	private:
		double originDepth_, radius_;
	};

}

#endif
