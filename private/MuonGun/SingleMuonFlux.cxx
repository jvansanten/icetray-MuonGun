
#include <MuonGun/SingleMuonFlux.h>

namespace I3MuonGun {

SingleMuonFlux::~SingleMuonFlux() {}

AdHocSingleMuonFlux::AdHocSingleMuonFlux() :
    a0_(0.003275), a1_(2.285), b0_(-0.09434), b1_(-0.3713),
    c0_(2.132), c1_(-0.1838), d0_(-0.0168), d1_(0.02342)
{}

double
AdHocSingleMuonFlux::operator()(double h, double x) const
{
	return a0_*pow(h, -a1_)*(x + (c0_ + c1_*h)*x*x)
	    *exp((b0_ + b1_*h)*(1./x + (d0_ + d1_*h)/(x*x)));
}

BMSSSingleMuonFlux::BMSSSingleMuonFlux() :
    k0a_(7.2e-3), k0b_(-1.927), k1a_(-0.581), k1b_(0.034)
{}

double
BMSSSingleMuonFlux::operator()(double h, double x) const
{
	return k0a_*pow(h, k0b_)*x*exp((k1a_*h + k1b_)/x);
}

}