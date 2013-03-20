/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_GENERATOR_H_INCLUDED
#define I3MUONGUN_GENERATOR_H_INCLUDED

#include <list>
#include <vector>
#include <boost/make_shared.hpp>

#include <icetray/I3PointerTypedefs.h>
#include <icetray/I3FrameObject.h>

class I3Particle;
class I3RandomService;
template <typename T> class I3Tree;
typedef I3Tree<I3Particle> I3MCTree;

namespace I3MuonGun {

/**
 * The radial offset and energy of each muon in a bundle
 */
struct BundleEntry {
	BundleEntry(double r=0., double e=0.) : radius(r), energy(e) {}
	double radius, energy;
	bool operator<(const BundleEntry &other) const
	{
		return other.energy < this->energy;
	}
	bool operator==(const BundleEntry &other) const
	{
		return (this->radius == other.radius) && (this->energy == other.energy);
	}
};
typedef std::list<BundleEntry> BundleConfiguration;

I3_FORWARD_DECLARATION(Surface);
I3_FORWARD_DECLARATION(SamplingSurface);
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
	
	void SetTotalEvents(double n) { numEvents_ = n; }
	double GetTotalEvents() const { return numEvents_; }
	
	/**
	 * @brief Calculate the differential number of events that
	 *        should have been generated by the represented scheme
	 *
	 * @param[in] axis   the bundle axis
	 * @param[in] bundle the radial offset and energy of each muon
	 *                   in the bundle
	 */
	double GetLogGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const;
	double GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const;
	
public:
	/**
	 * @brief Propose an injection surface for the given bundle configuration.
	 * 
	 * If variable-surface and fixed-surface generation schemes are to be combined,
	 * the variable surface must satisfy a few basic requirements:
	 * - The variable surface must always be *inside* any fixed surfaces
	 * - As the bundle properties are sampled at later and later times, the
	 *   proposed surface must retreat. For example, the surface must become
	 *   monotonically smaller for lower energies and multiplicities.
	 *
	 * @param[in] axis   an I3Particle representing the shower axis
	 * @param[in] bundle the radial offset and energy of each muon
	 *                   in the bundle
	 * @returns the surface where a bundle of the given characteristics
	 *          would be placed
	 */
	virtual SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const = 0;
	 
	/** Copy self into a shared pointer */
	virtual GenerationProbabilityPtr Clone() const = 0;
	
	/** @brief Compare to another GenerationProbability
	 *
	 * @returns true if the argument is identical to *this to
	 *          within a scale factor, false otherwise.
	 */
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const = 0;
protected:
	/**
	 * @brief Calculate the differential probability per event that the
	 *        given configuration was generated.
	 *
	 * For single muons, this is @f$ \log(dP/dE) [\log(1/GeV)]@f$, for bundles
	 * @f$ \log(d^2P/dEdr) [\log(1/GeV m)]@f$
	 *
	 * @param[in] axis   the bundle axis
	 * @param[in] bundle the radial offset and energy of each muon
	 *                   in the bundle
	 */
	virtual double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const = 0;

private:
	/** @brief The total number of events that should be generated */
	double numEvents_;
};

I3_POINTER_TYPEDEFS(GenerationProbability);

/**
 * @brief A collection of muon bundle generation schemes
 */
class GenerationProbabilityCollection : public GenerationProbability, public std::vector<GenerationProbabilityPtr> {
public:
	GenerationProbabilityCollection(GenerationProbabilityPtr, GenerationProbabilityPtr);
	void push_back(const GenerationProbabilityPtr&);
public:
	// GenerationProbability interface
	GenerationProbabilityPtr Clone() const;
	SamplingSurfaceConstPtr GetInjectionSurface(const I3Particle &axis, const BundleConfiguration &bundle) const;
	virtual bool IsCompatible(GenerationProbabilityConstPtr) const;
protected:
	/**
	 * Calculate the *total* probability that the given configuration was generated
	 * by any of the distributions in the colleciton.
	 */
	double GetLogGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
};

/** Scale the distribution by the given number of events */
GenerationProbabilityPtr operator*(double, GenerationProbabilityPtr);
GenerationProbabilityPtr operator*(GenerationProbabilityPtr, double);
GenerationProbabilityPtr operator*=(GenerationProbabilityPtr, double);
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
class Generator : public GenerationProbability, public I3FrameObject {
public:
	virtual ~Generator();
	/**
	 * @brief Generate a muon bundle.
	 *
	 * @param[in]  rng    A random number generator
	 * @param[out] tree   An I3MCTree to fill the generated bundle in to.
	 *                    The bundle axis should be used as the primary,
	 *                    with its type set to I3Particle::unknown
	 * @param[out] bundle the radial offset and energy of each muon
	 *                    in the bundle
	 */
	virtual void Generate(I3RandomService &rng, I3MCTree &tree, BundleConfiguration &bundle) const = 0;

	/**
	 * @brief Place a muon at a given radial offset and rotation with respect
	 * to the shower axis.
	 *
	 * @param[in] radius  perpendicular distance from the shower axis
	 * @param[in] azimuth rotation about the shower axis
	 * @param[in] surface place the muon on this surface, with timing
	 *                    adjusted so that it remains in the shower plane
	 * @param[in] axis    the shower axis
	 */
	static I3Particle CreateParallelTrack(double radius, double azimuth,
	    const Surface &surface, const I3Particle &axis);
};

I3_POINTER_TYPEDEFS(Generator);

}

#endif // I3MUONGUN_GENERATOR_H_INCLUDED
