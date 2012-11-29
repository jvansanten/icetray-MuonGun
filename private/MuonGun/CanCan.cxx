
#include <MuonGun/CanCan.h>
#include <MuonGun/I3MuonGun.h>
#include <dataclasses/I3Constants.h>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

namespace I3MuonGun {

CanCan::CanCan(double radius, double height, unsigned maxMultiplicity) :
    surface_(height, radius),
    originDepth_(I3Constants::SurfaceElev-I3Constants::OriginElev),
    maxMultiplicity_(maxMultiplicity), totalRate_(NAN)
{
	flux_ = boost::make_shared<AdHocSingleMuonFlux>();
	multiFraction_ = boost::make_shared<MultiplicityFraction>();
	radialDistribution_ = boost::make_shared<BMSSRadialDistribution>();
	energyDistribution_ = boost::make_shared<PotemkinEnergyDistribution>();
	
	// Any evaluated flux is guaranteed to be <= the maximum flux * the maximum differential area
	double thetaMax = atan(2*height_/(M_PI*radius_));
	maxFlux_ = (*flux_)(GetDepth(height_/2.), 1.)*GetDifferentialArea(cos(thetaMax));
	printf("Flux(%.1f,0) = %.3e\n", GetDepth(height_/2.), (*flux_)(GetDepth(height_/2.), 1.));
}

void
CanCan::SetRandomService(I3RandomServicePtr r)
{
	Distribution::SetRandomService(r);
	flux_->SetRandomService(r);
	multiFraction_->SetRandomService(r);
	radialDistribution_->SetRandomService(r);
	energyDistribution_->SetRandomService(r);
}

double
CanCan::GetDepth(double z) const
{
	return (originDepth_ - z)/I3Units::km;
}

double
CanCan::GetDifferentialArea(double coszen) const
{
	return M_PI*radius_*(radius_*coszen + (2*height_/M_PI)*sqrt(1-coszen*coszen));
}

// dAd Omega/dcos(theta) dphi (only one depth)
double
CanCan::GetDifferentialTopArea(double coszen) const
{
	return M_PI*radius_*(radius_*coszen);
}

// dAd Omega/dcos(theta) dphi dz (differential also in depth)
double
CanCan::GetDifferentialSideArea(double coszen) const
{
	return M_PI*radius_*((2/M_PI)*sqrt(1-coszen*coszen));
}

I3Position
CanCan::GetImpactPoint(const I3Direction &dir) const
{
	// The projection of a cylinder onto a plane whose
	// normal is inclined by `zenith` w.r.t to the cylinder
	// axis is a rectangle of width 2*r and height h*sin(theta)
	// capped with two half-ellipses of major axis r and
	// minor axis r*cos(theta). Pick a point from a uniform
	// distribution over this area.
	double a = sin(dir.GetZenith())*height_/2.;
	double b = cos(dir.GetZenith())*radius_;
	double x, y;
	do {
		x = radius_*rng_->Uniform(-1, 1);
		y = (a + b)*rng_->Uniform(-1, 1);
	} while (fabs(y) > a + b*sqrt(1 - (x*x)/(radius_*radius_)));
	I3Position impact = I3Position(y, x, 0);
	impact.RotateY(dir.GetZenith());
	impact.RotateZ(dir.GetAzimuth());
	
	// Now, project back to the entry point
	double l = surface_.GetIntersection(impact, dir).first;
	// std::cout << l << std::endl;
	impact.SetX(impact.GetX() + l*dir.GetX());
	impact.SetY(impact.GetY() + l*dir.GetY());
	impact.SetZ(impact.GetZ() + l*dir.GetZ());
	
	return impact;
}

std::pair<I3Particle, unsigned>
CanCan::GenerateAxis() const
{
	I3Particle p;
	unsigned multiplicity;
	double flux;
	do {
		// Pick a direction uniformly on the upper hemisphere
		p.SetDir(acos(rng_->Uniform()), rng_->Uniform(0., 2*M_PI));
		// Pick an impact point uniformly distributed in the projected area of the cylinder
		p.SetPos(GetImpactPoint(p.GetDir()));
		// Pick a number of muons in the bundle
		multiplicity = rng_->Integer(maxMultiplicity_)+1u;
		// Now, calculate the flux expectation at the chosen zenith angle
		// and at the depth where the shower axis crosses the cylinder
		double h = GetDepth(p.GetPos().GetZ());
		double zenith = p.GetDir().GetZenith();
		flux = (*flux_)(h, cos(zenith))
		    * (*multiFraction_)(h, cos(zenith), multiplicity)
		    * GetDifferentialArea(cos(zenith));
	} while (flux <= rng_->Uniform(0., maxFlux_));
	
	p.SetShape(I3Particle::Primary);
	p.SetLocationType(I3Particle::Anywhere);
	p.SetType(I3Particle::unknown);
	
	return std::make_pair(p, multiplicity);
}

double
CanCan::GenerateBundle(std::vector<I3Particle> &tracks) const
{
	// Pick a shower axis and number of muons to generate
	std::pair<I3Particle, unsigned> axis = GenerateAxis();
	tracks.push_back(axis.first);
	double depth = GetDepth(axis.first.GetPos().GetZ());
	double cos_theta = cos(axis.first.GetDir().GetZenith());
	
	double gen_prob = 2*M_PI*(*flux_)(depth, cos_theta)
	    * (*multiFraction_)(depth, cos_theta, axis.second)/GetTotalRate();
	
	// Draw N tracks from the radius/energy distribution at that depth and angle
	for (unsigned i=0; i < axis.second; i++) {
		Sample radius = radialDistribution_->Generate(depth, cos_theta, axis.second);
		Sample energy = energyDistribution_->Generate(depth, cos_theta, axis.second, radius.value);
		double azi = rng_->Uniform(0, 2*M_PI);
		I3Position offset(radius.value, 0, 0);
		offset.RotateY(axis.first.GetDir().GetZenith());
		offset.RotateZ(azi);
		
		I3Particle p;
		p.SetDir(axis.first.GetDir());
		p.SetPos(offset.GetX()+axis.first.GetPos().GetX(),
		         offset.GetY()+axis.first.GetPos().GetY(),
		         offset.GetZ()+axis.first.GetPos().GetZ());
		// TODO: shift the times and positions of off-axis tracks
		// so that they correspond to a plane wave crossing the sampling
		// surface at time 0
		p.SetEnergy(energy.value);
		p.SetType(I3Particle::MuMinus);
		p.SetShape(I3Particle::InfiniteTrack);
		tracks.push_back(p);
		
		gen_prob *= radius.prob*energy.prob;
	}
	
	return gen_prob;
}


namespace {
	
template <size_t N> struct multiply;

template <>
struct multiply<1> {
	typedef boost::function<double (double)> func_t;
	multiply(func_t f, func_t g, func_t h) : f_(f), g_(g), h_(h) {};
	func_t f_, g_, h_;
	double operator()(double x) const { return f_(x)*g_(x)*h_(x); }
};

template<>
struct multiply<2> {
	typedef boost::function<double (double, double)> func_t;
	multiply(func_t f, func_t g, func_t h) : f_(f), g_(g), h_(h) {};
	func_t f_, g_, h_;
	double operator()(double x, double y) const { return f_(x,y)*g_(x,y)*h_(x,y); }
};

}

// Calculate the expected rate dN/dt from the size of the sampling surface and
// the flux parametrization. The livetime is then N_generated/(dN/dt)
double
CanCan::CalculateTotalRate(unsigned multiplicity) const
{
	typedef boost::function<double (double)> f1;
	typedef boost::function<double (double, double)> f2;
	
	// First, integrate to find dN/dt on the top cap
	double top = 0;
	{
		f1 dN = boost::bind<double>(boost::cref(*flux_),
		    GetDepth(height_/2.), _1);
		f1 M = boost::bind<double>(boost::cref(*multiFraction_),
		    GetDepth(height_/2.), _1, multiplicity);
		f1 dOmega = boost::bind(&CanCan::GetDifferentialTopArea, this, _1);
		f1 dN_dOmega = multiply<1>(M, dN, dOmega);
		top = 2*M_PI*Integrate(dN_dOmega, 0., 1.);
	}
	
	// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
	double sides = 0.;
	{
		f2 dN = boost::bind<double>(boost::cref(*flux_),
		    boost::bind(&CanCan::GetDepth, this, _1), _2);
		f2 M = boost::bind<double>(boost::cref(*multiFraction_),
		    boost::bind(&CanCan::GetDepth, this, _1), _2, multiplicity);
		f2 dOmega = boost::bind(&CanCan::GetDifferentialSideArea, this, _2);
		f2 dN_dOmega = multiply<2>(M, dN, dOmega);
		double low[2] = {-height_/2., 0.};
		double high[2] = {height_/2., 1.};
		
		sides = 2*M_PI*Integrate(dN_dOmega, low, high, 2e-8, 2e-8, 10000u);
	}
	
	return top + sides;
}

double
CanCan::GetTotalRate() const
{
	if (std::isnan(totalRate_)) {
		totalRate_ = 0;
		for (unsigned m = 1; m <= maxMultiplicity_; m++)
			totalRate_ += CalculateTotalRate(m);
	}
	return totalRate_;
}

SingleMuonFlux::~SingleMuonFlux() {}

AdHocSingleMuonFlux::AdHocSingleMuonFlux() :
    a0_(0.001055), a1_(2.257), b0_(-0.08405), b1_(-0.4036),
    c0_(5.821), c1_(-0.5493), d0_(-0.01342), d1_(0.02279)
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

double
MultiplicityFraction::operator()(double depth, double cos_theta,
    unsigned multiplicity) const
{
	return 1.;
}

double
PotemkinEnergyDistribution::operator()(double depth, double cos_theta, 
    unsigned multiplicity, double radius, double energy) const
{
	return 1.;
}

Sample
PotemkinEnergyDistribution::Generate(double depth, double cos_theta, 
    unsigned multiplicity, double radius) const
{
	return Sample(1., 1.);
}

}