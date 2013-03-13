/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Generator.h>
#include <MuonGun/Surface.h>

#include <icetray/I3Module.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/I3Double.h>
#include <phys-services/I3RandomService.h>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


namespace I3MuonGun {

GenerationProbability::~GenerationProbability() {}
Generator::~Generator() {}

double
GenerationProbability::GetLogGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return std::log(double(numEvents_)) + GetLogGenerationProbability(axis, bundle);
}

double
GenerationProbability::GetGeneratedEvents(const I3Particle &axis, const BundleConfiguration &bundle) const
{
	return numEvents_*std::exp(GetLogGenerationProbability(axis, bundle));
}

double
GenerationProbability::GetLogGeneratedBundles(const I3Particle &axis, unsigned multiplicity) const
{
	return std::log(double(numEvents_)) + GetLogBundleGenerationProbability(axis, multiplicity);
}

GenerationProbabilityCollection::GenerationProbabilityCollection(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{
	push_back(p1);
	push_back(p2);
}

double
GenerationProbabilityCollection::GetLogGenerationProbability(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	// Collect log probabilities from members
	std::vector<double> values(this->size(), -std::numeric_limits<double>::infinity());
	BOOST_FOREACH(const value_type &p, *this)
		if (p)
			values.push_back(p->GetLogGeneratedEvents(axis, bundle));
	
	// Calculate log(sum(exp)) in a numerically stable way
	double bias = *std::max_element(values.begin(), values.end());
	double prob = 0.;
	BOOST_FOREACH(double v, values)
		prob += std::exp(v-bias);
	
	return bias + std::log(prob);
}

/**
 * @brief Find the inner-most injection surface in the collection
 *
 * Generators are allowed to use different sampling surfaces, and even
 * to change the size and shape of their injection surfaces based on the
 * properties of the muon bundle. The combined weighting for simulations
 * on different injection surfaces only makes sense, however, if the flux
 * is measured on the inner-most surface. For a given shower axis, this
 * is the one whose entry point is furthest away.
 */
SamplingSurfaceConstPtr
GenerationProbabilityCollection::GetInjectionSurface(const I3Particle &axis,
    const BundleConfiguration &bundle) const
{
	std::map<double, SamplingSurfaceConstPtr> surfaces;
	BOOST_FOREACH(const value_type &p, *this) {
		SamplingSurfaceConstPtr surface = p->GetInjectionSurface(axis, bundle);
		std::pair<double, double> steps =
		    surface->GetIntersection(axis.GetPos(), axis.GetDir());
		// If any surface is missed the entire event should be
		// weighted to 0, so we ensure that any that are missed
		// will be considered "innermost"
		surfaces[std::isfinite(steps.first) ? steps.first : std::numeric_limits<double>::infinity()] = surface;
	}
	
	if (surfaces.size() == 0)
		log_fatal("Empty collection!");
	return (--surfaces.end())->second;
}

GenerationProbabilityPtr
GenerationProbabilityCollection::Clone() const
{
	return boost::make_shared<GenerationProbabilityCollection>(*this);
}

GenerationProbabilityPtr
operator*=(GenerationProbabilityPtr p, size_t n)
{
	p->SetTotalEvents(p->GetTotalEvents()*n);
	
	return p;
}

GenerationProbabilityPtr
operator*(GenerationProbabilityPtr p, size_t n)
{
	GenerationProbabilityPtr pn = p->Clone();
	pn *= n;
	
	return pn;
}

GenerationProbabilityPtr
operator*(size_t n, GenerationProbabilityPtr p)
{
	return p*n;
}

GenerationProbabilityPtr
operator+(GenerationProbabilityPtr p1, GenerationProbabilityPtr p2)
{
	boost::shared_ptr<GenerationProbabilityCollection> c1 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p1);
	boost::shared_ptr<GenerationProbabilityCollection> c2 =
	    boost::dynamic_pointer_cast<GenerationProbabilityCollection>(p2);
	if (c1 && c2) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		std::copy(c2->begin(), c2->end(), std::back_inserter(*c1));
		return c1;
	} else if (c1) {
		c1 = boost::make_shared<GenerationProbabilityCollection>(*c1);
		c1->push_back(p2);
		return c1;
	} else if (c2) {
		c2 = boost::make_shared<GenerationProbabilityCollection>(*c2);
		c2->push_back(c1);
		return c2;
	}
	
	return boost::make_shared<GenerationProbabilityCollection>(p1, p2);
}

namespace {

inline I3Position
operator+(const I3Position &a, const I3Position &b)
{
	return I3Position(a.GetX()+b.GetX(), a.GetY()+b.GetY(), a.GetZ()+b.GetZ());
}

inline I3Position&
operator+=(I3Position &a, const I3Position &b)
{
	a.SetX(a.GetX() + b.GetX());
	a.SetY(a.GetY() + b.GetY());
	a.SetZ(a.GetZ() + b.GetZ());
	return a;
}

inline I3Position
operator*(double v, const I3Direction &d)
{
	return I3Position(v*d.GetX(), v*d.GetY(), v*d.GetZ());
}

}

I3Particle
Generator::CreateParallelTrack(double radius, double azimuth,
    const Surface &surface, const I3Particle &axis)
{
	I3Particle track;
	track.SetLocationType(I3Particle::InIce);
	track.SetType(I3Particle::MuMinus);
	track.SetDir(axis.GetDir());
	track.SetSpeed(I3Constants::c);
	track.SetPos(axis.GetPos());
	track.SetTime(axis.GetTime());
	
	if (radius > 0) {
		// Shift the track parallel to the axis
		I3Position offset(radius, 0, 0);
		offset.RotateY(axis.GetDir().GetZenith());
		offset.RotateZ(azimuth);
		offset += axis.GetPos();
		// Find the distance from the offset position to the sampling
		// surface, and shift so that all tracks have their origin
		// on the surface, but remain in a plane
		double shift = 
		    surface.GetIntersection(offset, track.GetDir()).first
		    -surface.GetIntersection(axis.GetPos(), axis.GetDir()).first;
		if (std::isfinite(shift)) {
			track.SetTime(axis.GetTime() + shift/track.GetSpeed());
			track.SetPos(offset + shift*track.GetDir());
		} else {
			track.SetPos(offset);
		}
	}
	
	return track;
}

/**
 * @brief Interface between Generator and IceTray
 */
class GeneratorModule : public I3Module {
public:
	GeneratorModule(const I3Context &ctx) : I3Module(ctx), maxEvents_(0), numEvents_(0)
	{
		AddOutBox("OutBox");
		AddParameter("Generator", "Muon bundle generator", generator_);
		
		mctreeName_ = "I3MCTree";
	}
	
	void Configure()
	{
		GetParameter("Generator", generator_);
		
		rng_ = context_.Get<I3RandomServicePtr>();
		if (!rng_)
			log_fatal("No RandomService configured!");
		maxEvents_ = generator_->GetTotalEvents();
		
		firstFrame_ = true;
	}
	
	void DAQ(I3FramePtr frame)
	{
		if (firstFrame_) {
			firstFrame_ = false;
			I3FramePtr sframe = boost::make_shared<I3Frame>('S');
			sframe->Put(GetName(), generator_);
			PushFrame(sframe);
		}
		
		I3MCTreePtr mctree = boost::make_shared<I3MCTree>();
		BundleConfiguration bundlespec;
		
		generator_->Generate(*rng_, *mctree, bundlespec);
		
		frame->Put(mctreeName_, mctree);
		
		PushFrame(frame);
		if (++numEvents_ >= maxEvents_)
			RequestSuspension();
	}
private:
	GeneratorPtr generator_;
	I3RandomServicePtr rng_;
	size_t maxEvents_, numEvents_;
	std::string mctreeName_;
	bool firstFrame_;
};

}

I3_MODULE(I3MuonGun::GeneratorModule);
