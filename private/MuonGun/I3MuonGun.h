
#include <boost/function.hpp>

class I3Position;
class I3Direction;

// Commonly-used bits
namespace I3MuonGun {

// Find the distances to the points of intersection with a centered at (0,0,0) and aligned along the z axis.
std::pair<double, double> CylinderIntersection(const I3Position &p, const I3Direction &dir, double length, double radius);

double Integrate(boost::function<double (double)> f, double low, double high, double epsabs=1.49e-8, double epsrel=1.49e-8, size_t limit=50);
double Integrate(boost::function<double (double, double)> f, double low[2], double high[2], double epsabs=1.49e-8, double epsrel=1.49e-8, size_t limit=0);


}