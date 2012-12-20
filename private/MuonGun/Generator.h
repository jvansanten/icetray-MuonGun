
#ifndef I3MUONGUN_GENERATOR_H_INCLUDED
#define I3MUONGUN_GENERATOR_H_INCLUDED

#include <MuonGun/Distribution.h>
#include <boost/make_shared.hpp>

class I3Particle;
template <typename T> class I3Tree;
typedef I3Tree<I3Particle> I3MCTree;

namespace I3MuonGun {

typedef std::vector<std::pair<double, double> > BundleConfiguration;

I3_FORWARD_DECLARATION(Surface);
I3_FORWARD_DECLARATION(GenerationProbability);

class GenerationProbability {
public:
	GenerationProbability() : numEvents_(1) {}
	virtual ~GenerationProbability();
	
	void SetTotalEvents(size_t n) { numEvents_ = n; }
	size_t GetTotalEvents() const { return numEvents_; }
	
	double GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const;
public:
	/// Copy self into a shared pointer
	virtual GenerationProbabilityPtr Clone() const = 0;
protected:
	/// Calculate the probability that the given configuration was generated
	virtual double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const = 0;

private:
	size_t numEvents_;
};

I3_POINTER_TYPEDEFS(GenerationProbability);

class GenerationProbabilityCollection : public GenerationProbability, public std::vector<GenerationProbabilityPtr> {
public:
	GenerationProbabilityCollection(GenerationProbabilityPtr, GenerationProbabilityPtr);
	GenerationProbabilityPtr Clone() const;
protected:
	double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const;
};

GenerationProbabilityPtr operator*(size_t, GenerationProbabilityPtr);
GenerationProbabilityPtr operator*(GenerationProbabilityPtr, size_t);
GenerationProbabilityPtr operator*=(GenerationProbabilityPtr, size_t);
GenerationProbabilityPtr operator+(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2);

class Generator : public Distribution, public GenerationProbability {
public:
	virtual ~Generator();
	virtual void Generate(I3MCTree &tree, BundleConfiguration &bundle) const = 0;
	virtual SurfaceConstPtr GetInjectionSurface() const = 0;
};

I3_POINTER_TYPEDEFS(Generator);

}

#endif // I3MUONGUN_GENERATOR_H_INCLUDED
