
#include <MuonGun/CanCan.h>
#include <MuonGun/I3MuonGun.h>
#include <dataclasses/I3Constants.h>

#include <boost/bind.hpp>

namespace I3MuonGun {

CanCan::CanCan(double radius, double height) :
    flux_(), density_(0.917), radius_(radius), height_(height),
    originDepth_(I3Constants::SurfaceElev-I3Constants::OriginElev) 
{
	double thetaMax = atan(2*height_/(M_PI*radius_));
	// Any evaluated flux is guaranteed to be <= the maximum flux * the maximum differential area
	maxFlux_ = flux_.GetFlux(GetDepth(height_/2.), 1.)*GetDifferentialArea(cos(thetaMax));
	printf("Flux(%.1f,0) = %.3e\n", GetDepth(height_/2.), flux_.GetFlux(GetDepth(height_/2.), 1.));
}

double
CanCan::GetDepth(double z) const
{
	// Top 200 m are less dense
	return (200*density_*0.832 + (originDepth_ - 200 - z)*density_*1.005)/I3Units::km;
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
	double l = CylinderIntersection(impact, dir, height_, radius_).first;
	// std::cout << l << std::endl;
	impact.SetX(impact.GetX() + l*dir.GetX());
	impact.SetY(impact.GetY() + l*dir.GetY());
	impact.SetZ(impact.GetZ() + l*dir.GetZ());
	
	return impact;
}

I3Particle CanCan::Generate() const
{
	I3Particle p;
	double flux;
	do {
		// Pick a direction uniformly on the upper hemisphere
		p.SetDir(acos(rng_->Uniform()), rng_->Uniform(0., 2*M_PI));
		// Pick an impact point uniformly distributed in the projected area of the cylinder
		p.SetPos(GetImpactPoint(p.GetDir()));
		// Now, calculate the flux expectation at the chosen zenith angle
		// and at the depth where the track crosses the cylinder
		double h = GetDepth(p.GetPos().GetZ());
		double zenith = p.GetDir().GetZenith();
		flux = flux_.GetFlux(h, cos(zenith))*GetDifferentialArea(cos(zenith));
	} while (flux <= rng_->Uniform(0., maxFlux_));
	
	return p;
}

namespace {
	
template <size_t N> struct multiply;

template <>
struct multiply<1> {
	typedef boost::function<double (double)> func_t;
	multiply(func_t f, func_t g) : f_(f), g_(g) {};
	func_t f_, g_;
	double operator()(double x) const { return f_(x)*g_(x); }
};

template<>
struct multiply<2> {
	typedef boost::function<double (double, double)> func_t;
	multiply(func_t f, func_t g) : f_(f), g_(g) {};
	func_t f_, g_;
	double operator()(double x, double y) const { return f_(x,y)*g_(x,y); }
};

}

// Calculate the expected rate dN/dt from the size of the sampling surface and
// the flux parametrization. The livetime is then N_generated/(dN/dt)
double
CanCan::GetLivetime() const
{
	typedef boost::function<double (double)> f1;
	typedef boost::function<double (double, double)> f2;
	
	// First, integrate to find dN/dt on the top cap
	double top = 0;
	{
		f1 dN = boost::bind(&SingleMuonFlux::GetFlux, &flux_, GetDepth(height_/2.), _1);
		f1 dOmega = boost::bind(&CanCan::GetDifferentialTopArea, this, _1);
		f1 dN_dOmega = multiply<1>(dN, dOmega);
		top = 2*M_PI*Integrate(dN_dOmega, 0., 1.);
	}
	
	// Now, the more complicated bit: integrate over the sides. The flux is now a function of both depth and zenith!
	double sides = 0.;
	{
		f2 dN = boost::bind(&SingleMuonFlux::GetFlux, &flux_, boost::bind(&CanCan::GetDepth, this, _1), _2);
		f2 dOmega = boost::bind(&CanCan::GetDifferentialSideArea, this, _2);
		f2 dN_dOmega = multiply<2>(dN, dOmega);
		double low[2] = {-height_/2., 0.};
		double high[2] = {height_/2., 1.};
		
		sides = 2*M_PI*Integrate(dN_dOmega, low, high, 2e-8, 2e-8, 10000u);
	}
	
	return top + sides;
}

SingleMuonFlux::SingleMuonFlux() :
    a0_(0.001731), a1_(2.202), b0_(-0.0808), b1_(-0.3866),
    c0_(2.229), c1_(-0.2239), d0_(-0.01059), d1_(0.02661)
{}

double
SingleMuonFlux::GetFlux(double h, double x) const
{
	return a0_*pow(h, -a1_)*(x + (c0_ + c1_*h)*x*x)
	    *exp((b0_ + b1_*h)*(1./x + (d0_ + d1_*h)/(x*x)));
}

}