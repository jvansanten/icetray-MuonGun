/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_FLUX_H_INCLUDED
#define I3MUONGUN_FLUX_H_INCLUDED

#include <photospline/I3SplineTable.h>

namespace I3MuonGun {

/**
 * @brief A parameterization of the total muon-bundle flux
 * 
 * The total flux (in units of @f$ [ m^{-2} sr^{-2} s^{-1} ]@f$)
 * is parameterized in terms of vertical depth @f$ km @f$, the
 * cosine of the zenith angle, and bundle multiplicity.
 */
class Flux {
public:
	Flux();
	virtual ~Flux();
	typedef double result_type;
	// TODO: short-circuit if multiplicity outside bounds
	virtual double operator()(double depth, double cos_theta, unsigned multiplicity) const = 0;

	unsigned GetMaxMultiplicity() const { return maxMultiplicity_; }
	unsigned GetMinMultiplicity() const { return minMultiplicity_; }
	
	void SetMaxMultiplicity(unsigned m) { maxMultiplicity_ = m; }
	void SetMinMultiplicity(unsigned m) { minMultiplicity_ = m; }

private:
	unsigned minMultiplicity_, maxMultiplicity_;
};

I3_POINTER_TYPEDEFS(Flux);

/**
 * @brief Total flux according to Becherini et al.
 */
class BMSSFlux : public Flux {
public:
	BMSSFlux();
	double operator()(double depth, double cos_theta, unsigned multiplicity) const;
private:
	double k0a_, k0b_, k1a_, k1b_;
	double v0a_, v0b_, v0c_, v1a_, v1b_;
};

/**
 * @brief Total flux, fit to a tensor-product B-spline surface
 */
class SplineFlux : public Flux {
public:
	SplineFlux(const std::string &singles, const std::string &bundles);
	double operator()(double depth, double cos_theta, unsigned multiplicity) const;
private:
	I3SplineTable singles_;
	I3SplineTable bundles_;
};

}

#endif // I3MUONGUN_FLUX_H_INCLUDED