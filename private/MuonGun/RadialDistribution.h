/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_RADIALDISTRIBUTION_H
#define I3MUONGUN_RADIALDISTRIBUTION_H

#include <photospline/I3SplineTable.h>
#include <icetray/I3PointerTypedefs.h>

class I3Position;
class I3RandomService;

namespace I3MuonGun {

/**
 * @brief The distribution of distance between a muon and the bundle axis
 */
class RadialDistribution {
public:
	virtual ~RadialDistribution();
	// TODO: short-circuit for multiplicity == 1
	/**
	 * Calculate the probability of obtaining the given radial offset
	 *
	 * @param[in] depth        vertical depth in km
	 * @param[in] cos_theta    cosine of zenith angle
	 * @param[in] multiplicity number of muons in the bundle
	 * @param[in] radius       distance to bundle axis
	 * @returns a properly normalized propability @f$ dP/dr \,\, [m^{-1}] @f$
	 */
	virtual double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const = 0;
	/**
	 * Draw a sample from the distribution of radii
	 *
	 * @param[in] depth        vertical depth in km
	 * @param[in] cos_theta    cosine of zenith angle
	 * @param[in] multiplicity number of muons in the bundle
	 * @returns a radius in m
	 */
	virtual double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const = 0;
};

I3_POINTER_TYPEDEFS(RadialDistribution);

/**
 * @brief Radial distribution according to Becherini et al.
 */
class BMSSRadialDistribution : public RadialDistribution {
public:
	BMSSRadialDistribution();
	double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const;
private:
	double GetMeanRadius(double, double, unsigned) const;
	double GetShapeParameter(double, double, unsigned) const;
	double GetGenerationProbability(double, double, double) const;
	double rho0a_, rho0b_, rho1_, theta0_, f_, alpha0a_, alpha0b_, alpha1a_, alpha1b_;
	double rmax_;
};

I3_POINTER_TYPEDEFS(BMSSRadialDistribution);

/**
 * @brief Radial distribution fit to a tensor-product B-spline surface
 *
 * The surface is fit to @f$ d \log{P} / d{r^2} @f$ to remove the factor
 * of differential area implicit in @f$ dP / dr @f$
 */
class SplineRadialDistribution : public RadialDistribution, private I3SplineTable {
public:
	SplineRadialDistribution(const std::string&);
	double operator()(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
	double Generate(I3RandomService &rng, double depth, double cos_theta,
	    unsigned multiplicity) const;
private:
	double RawProbability(double depth, double cos_theta,
	    unsigned multiplicity, double radius) const;
};

}

#endif // I3MUONGUN_RADIALDISTRIBUTION_H
