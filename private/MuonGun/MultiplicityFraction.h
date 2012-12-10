
#ifndef I3MUONGUN_MULTIPLICITYFRACTION_H_INCLUDED
#define I3MUONGUN_MULTIPLICITYFRACTION_H_INCLUDED

#include <MuonGun/Distribution.h>

namespace I3MuonGun {
	
class MultiplicityFraction : public Distribution {
public:
	MultiplicityFraction();
	virtual ~MultiplicityFraction();
	typedef double result_type;
	virtual double operator()(double depth, double cos_theta, unsigned multiplicity) const = 0;
	
	unsigned GetMax() const { return max_; }
	unsigned GetMin() const { return min_; }
	void SetMax(unsigned v) { max_ = v; }
	void SetMin(unsigned v) { min_ = v; }
private:
	unsigned min_, max_;
};

I3_POINTER_TYPEDEFS(MultiplicityFraction);

class BMSSMultiplicityFraction : public MultiplicityFraction {
public:
	BMSSMultiplicityFraction();
	virtual double operator()(double depth, double cos_theta, unsigned multiplicity) const;
private:
	double v0a_, v0b_, v0c_, v1a_, v1b_;
};

}

#endif