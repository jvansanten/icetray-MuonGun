
#include <MuonGun/MultiplicityFraction.h>

namespace I3MuonGun {

MultiplicityFraction::MultiplicityFraction() : min_(1), max_(1) {}

MultiplicityFraction::~MultiplicityFraction() {};

BMSSMultiplicityFraction::BMSSMultiplicityFraction() :
    v0a_(0.01041), v0b_(0.09912), v0c_(2.712), v1a_(0.01615), v1b_(0.6010) {}

double
BMSSMultiplicityFraction::operator()(double h, double ct,
    unsigned m) const
{
	return pow(m, -(v0a_*h*h + v0b_*h + v0c_)*exp(v1a_*exp(v1b_*h)/ct));
}

}