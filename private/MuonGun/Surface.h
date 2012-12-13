
#ifndef I3MUONGUN_SURFACE_H_INCLUDED
#define I3MUONGUN_SURFACE_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>
#include <boost/function.hpp>

class I3Position;
class I3Direction;
class I3RandomService;

namespace I3MuonGun {

	class Surface {
	public:
		virtual ~Surface();
		virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const = 0;
	};
	
	I3_POINTER_TYPEDEFS(Surface);
	
	class SamplingSurface : public Surface {
	public:
		virtual ~SamplingSurface();
		// Get d(A_projected)/d(cos(theta)) at a particular cos(theta)
		virtual double GetDifferentialArea(double coszen) const = 0;
		// Get the maximum of d(A_projected)/d(cos(theta))
		virtual double GetMaxDifferentialArea() const = 0;
		virtual double GetMinDepth() const = 0;
	
		// Integrate a flux (defined in terms of dN/dOmega(depth [km], cos(theta)))
		// over the outer surface
		virtual double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const = 0;
	
		// Sample an impact point and direction from an isotropic flux between cosMin and cosMax,
		// and return the differential area in the chosen direction
		virtual double SampleImpactRay(I3Position &pos, I3Direction &dir, I3RandomService &rng, double cosMin=0, double cosMax=1) const = 0;
	};
	
	I3_POINTER_TYPEDEFS(SamplingSurface);

	class Cylinder : public SamplingSurface {
	public:
		Cylinder(double length, double radius) : length_(length), radius_(radius) {};
		std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	
		// Get d(A_projected)/d(cos(theta)) at a particular cos(theta)
		double GetDifferentialArea(double coszen) const;
		// Get the maximum of d(A_projected)/d(cos(theta))
		double GetMaxDifferentialArea() const;
		double GetMinDepth() const;
	
		// Integrate a flux (defined in terms of dN/dOmega(depth [km], cos(theta)))
		// over the outer surface
		double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const;
	
		// Sample an impact point and direction from an isotropic flux between cosMin and cosMax,
		// and return the differential area in the chosen direction
		double SampleImpactRay(I3Position &pos, I3Direction &dir, I3RandomService &rng, double cosMin=0, double cosMax=1) const;
	private:
		double GetDifferentialTopArea(double cos_zenith) const;
		double GetDifferentialSideArea(double cos_zenith) const;
	
		double length_, radius_;
	};

	I3_POINTER_TYPEDEFS(Cylinder);

	class Sphere : public Surface {
	public:
		Sphere(double originDepth, double radius) : originDepth_(originDepth), radius_(radius) {};
		std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	private:
		double originDepth_, radius_;
	};

}

#endif
