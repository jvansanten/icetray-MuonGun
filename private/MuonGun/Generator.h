/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * @version $Revision$
 * @date $Date$
 */

#ifndef I3MUONGUN_GENERATOR_H_INCLUDED
#define I3MUONGUN_GENERATOR_H_INCLUDED

#include <vector>
#include <boost/make_shared.hpp>

#include <icetray/I3PointerTypedefs.h>

class I3Particle;
class I3RandomService;
template <typename T> class I3Tree;
typedef I3Tree<I3Particle> I3MCTree;

namespace I3MuonGun {

/**
 * The radial offset and energy of each muon in a bundle
 */
typedef std::vector<std::pair<double, double> > BundleConfiguration;

I3_FORWARD_DECLARATION(Surface);
I3_FORWARD_DECLARATION(GenerationProbability);

/**
 * @brief A muon bundle generation scheme
 *
 * GenerationProbability represents the normalization required for WeightCalculator
 */
class GenerationProbability {
public:
	GenerationProbability() : numEvents_(1) {}
	virtual ~GenerationProbability();
	
	void SetTotalEvents(size_t n) { numEvents_ = n; }
	size_t GetTotalEvents() const { return numEvents_; }
	
	/**
	 * Calculate the number of events with the given configuration that
	 * that should have been produced
	 *
	 * @param[in] axis   an I3Particle representing the shower axis
	 * @param[in] bundle the radial offset and energy of each muon
	 *                   in the bundle
	 */
	double GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const;
public:
	/** Copy self into a shared pointer */
	virtual GenerationProbabilityPtr Clone() const = 0;
protected:
	/**
	 * Calculate the probability that the given configuration was generated
	 *
	 * @param[in] axis   an I3Particle representing the shower axis
	 * @param[in] bundle the radial offset and energy of each muon
	 *                   in the bundle
	 */
	virtual double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const = 0;

private:
	/** The total number of events that should be generated */
	size_t numEvents_;
};

I3_POINTER_TYPEDEFS(GenerationProbability);

/**
 * @brief A collection of muon bundle generation schemes
 */
class GenerationProbabilityCollection : public GenerationProbability, public std::vector<GenerationProbabilityPtr> {
public:
	GenerationProbabilityCollection(GenerationProbabilityPtr, GenerationProbabilityPtr);
	GenerationProbabilityPtr Clone() const;
protected:
	/**
	 * Calculate the *total* probability that the given configuration was generated
	 * by any of the distributions in the colleciton.
	 */
	double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
};

/** Scale the distribution by the given number of events */
GenerationProbabilityPtr operator*(size_t, GenerationProbabilityPtr);
GenerationProbabilityPtr operator*(GenerationProbabilityPtr, size_t);
GenerationProbabilityPtr operator*=(GenerationProbabilityPtr, size_t);
/**
 * Combine the distributions to form a GenerationProbabilityCollection
 * (or append to it if any of the arguments is already one)
 */
GenerationProbabilityPtr operator+(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2);

/**
 * @brief A muon bundle generator
 *
 * Generators know both how to draw muon bundles from some distribution and
 * how to calculate the probability that they would have drawn some arbitrary
 * bundle from that same distribution.
 */
class Generator : public GenerationProbability {
public:
	virtual ~Generator();
	/**
	 * Generate a muon bundle.
	 *
	 * @param[in]  rng    A random number generator
	 * @param[out] tree   An I3MCTree to fill the generated bundle in to.
	 *                    The bundle axis should be used as the primary,
	 *                    with its type set to I3Particle::unknown
	 * @param[out] bundle the radial offset and energy of each muon
	 *                    in the bundle
	 */
	virtual void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const = 0;
	/** Return the surface where bundles will be placed */
	virtual SurfaceConstPtr GetInjectionSurface() const = 0;
};

I3_POINTER_TYPEDEFS(Generator);

}

#endif // I3MUONGUN_GENERATOR_H_INCLUDED
