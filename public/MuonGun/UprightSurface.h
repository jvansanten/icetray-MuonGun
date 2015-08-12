
#ifndef I3MUONGUN_UPRIGHTSURFACE_H_INCLUDED
#define I3MUONGUN_UPRIGHTSURFACE_H_INCLUDED

#include <MuonGun/I3MuonGun.h>
#include <boost/bind.hpp>

namespace {

inline double
integrate_area(double a, double b, double cap, double sides)
{
	return 2*M_PI*(cap*(b*b-a*a) +
	    (sides/2.)*(acos(a) - acos(b) + sqrt(1-a*a)*a + sqrt(1-b*b)*b));
}

}

namespace I3MuonGun { namespace detail {
/**
 * @brief A surface consisting only of vertical and horizontal faces
 */
template <typename Base>
class UprightSurface : public Base {
public:
	// SamplingSurface interface
	double GetAcceptance(double cosMin=0, double cosMax=1) const
	{
		double cap = GetTopArea();
		double sides = GetSideArea();
		if (cosMin >= 0 && cosMax >= 0)
			return integrate_area(cosMin, cosMax, cap, sides);
		else if (cosMin < 0 && cosMax <= 0)
			return integrate_area(-cosMax, -cosMin, cap, sides);
		else if (cosMin < 0 && cosMax > 0)
			return integrate_area(0, -cosMin, cap, sides)
			    + integrate_area(0, cosMax, cap, sides);
		else
			log_fatal("Can't deal with zenith range [%.1e, %.1e]", cosMin, cosMax);
		return NAN;
	}
	double GetMinDepth() const
	{
		return GetDepth(GetZRange().second);
	}

	double IntegrateFlux(boost::function<double (double, double)> flux, double cosMin=0, double cosMax=1) const
	{
		typedef boost::function<double (double)> f1;
		typedef boost::function<double (double, double)> f2;
	
		double total = 0;
		std::pair<double, double> z_range = GetZRange();
	
		// First, integrate to find dN/dt on the cap(s)
		{
			f1 dN = boost::bind<double>(flux, GetDepth(z_range.second), _1);
			f1 dOmega = boost::bind(&UprightSurface<Base>::GetDifferentialTopArea, this, _1);
			f1 dN_dOmega = detail::multiply<1>(dN, dOmega);
			total += 2*M_PI*Integrate(dN_dOmega, cosMin, cosMax, 1e-3, 1e-3);
		}
	
		// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
		{
			f2 dN = boost::bind(flux, boost::bind(GetDepth, _1), _2);
			f2 dOmega = boost::bind(&UprightSurface<Base>::GetDifferentialSideArea, this, _2);
			f2 dN_dOmega = detail::multiply<2>(dN, dOmega);
			boost::array<double, 2> low = {{z_range.first, cosMin}};
			boost::array<double, 2> high = {{z_range.second, cosMax}};
			total += 2*M_PI*Integrate(dN_dOmega, low, high, 1e-3, 1e-3, 10000u);
		}
	
		return total;
	}

protected:
	virtual double GetSideArea() const = 0;
	virtual double GetTopArea() const = 0;
	virtual double GetLength() const = 0;
	virtual std::pair<double, double> GetZRange() const = 0;

	double GetDifferentialTopArea(double coszen) const
	{
		return std::abs(coszen)*GetTopArea();
	}
	double GetDifferentialSideArea(double coszen) const
	{
		return GetSideArea()/GetLength()*sqrt(1-coszen*coszen);
	}

protected:
	// not a generic way to forward ctors, but good enough for this use
	template <typename A, typename B>
	UprightSurface(A a, B b) : Base(a,b) {};
	template <typename A, typename B, typename C>
	UprightSurface(A a, B b, C c) : Base(a,b,c) {};
	
	UprightSurface() {};

private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &ar, unsigned version)
	{
		if (version > 0)
			log_fatal_stream("Version "<<version<<" is from the future");
	
		ar & make_nvp("Base", base_object<Base>(*this));
	}
};

}}

#endif
