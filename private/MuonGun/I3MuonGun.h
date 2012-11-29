
#ifndef MUONGUN_I3MUONGUN_H_INCLUDED
#define MUONGUN_I3MUONGUN_H_INCLUDED

#include <boost/function.hpp>

class I3Position;
class I3Direction;

// Commonly-used bits
namespace I3MuonGun {

class Surface {
public:
	virtual ~Surface();
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const = 0;
};

class Cylinder : public Surface {
public:
	Cylinder(double length, double radius) : length_(length), radius_(radius) {};
	std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
private:
	double length_, radius_;
};

class Sphere : public Surface {
public:
	Sphere(double originDepth, double radius) : originDepth_(originDepth), radius_(radius) {};
	std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
private:
	double originDepth_, radius_;
};

double Integrate(boost::function<double (double)> f, double low, double high, double epsabs=1.49e-8, double epsrel=1.49e-8, size_t limit=50);
double Integrate(boost::function<double (double, double)> f, double low[2], double high[2], double epsabs=1.49e-8, double epsrel=1.49e-8, size_t limit=0);


}

#endif
