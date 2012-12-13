
#ifndef I3MUONGUN_GENERATOR_H_INCLUDED
#define I3MUONGUN_GENERATOR_H_INCLUDED

#include <MuonGun/Distribution.h>
#include <dataclasses/physics/I3MCTree.h>

namespace I3MuonGun {

typedef std::vector<std::pair<double, double> > BundleConfiguration;

class Generator : public Distribution {
public:
	Generator() : numEvents_(1) {}
	virtual ~Generator();
	
	void SetTotalEvents(size_t n) { numEvents_ = n; }
	size_t GetTotalEvents() const { return numEvents_; }
	
	double GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const;
	
	virtual void Generate(I3MCTree &tree, BundleConfiguration &bundle) const = 0;
protected:
	virtual double GetGenerationProbability(const I3Particle &axis, const BundleConfiguration &bundle) const = 0;
private:
	size_t numEvents_;
};

I3_POINTER_TYPEDEFS(Generator);

}

#endif // I3MUONGUN_GENERATOR_H_INCLUDED
