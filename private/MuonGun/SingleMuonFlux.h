
#ifndef I3MUONGUN_SINGLEMUONFLUX_H_INCLUDED
#define I3MUONGUN_SINGLEMUONFLUX_H_INCLUDED

#include <MuonGun/Distribution.h>

namespace I3MuonGun {

class SingleMuonFlux : public Distribution {
public:
	virtual ~SingleMuonFlux();
	typedef double result_type;
	virtual double operator()(double depth, double cos_theta) const = 0;
};

I3_POINTER_TYPEDEFS(SingleMuonFlux);

// From Becherini et al.
class BMSSSingleMuonFlux : public SingleMuonFlux {
public:
	BMSSSingleMuonFlux();
	double operator()(double depth, double cos_theta) const;
private:
	double k0a_, k0b_, k1a_, k1b_;
};

// Becherini et al, re-fit to CORSIKA + MMC simulation
class AdHocSingleMuonFlux : public SingleMuonFlux {
public:
	AdHocSingleMuonFlux();
	double operator()(double depth, double cos_theta) const;
private:
	double a0_, a1_, b0_, b1_, c0_, c1_, d0_, d1_;
};

}

#endif
