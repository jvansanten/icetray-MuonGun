#include <I3Test.h>

#include <phys-services/I3Calculator.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3Units.h>

TEST_GROUP(DepthTest)

I3Direction
RotateToZenith(const I3Direction &direction, const I3Direction &dir)
{
	I3Direction p(dir);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Direction
RotateToZenith(const I3Direction &direction, const I3Position &pos)
{
	I3Direction p(pos);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Particle
RotateToZenith(const I3Direction &direction, I3Particle &part)
{
	I3Particle p(part);
	p.SetDir(RotateToZenith(direction, p.GetDir()));
	p.SetPos(RotateToZenith(direction, p.GetPos()));
	return p;
}

double
Angle(const I3Direction &d1, const I3Direction &d2)
{
	return acos(d1.GetX()*d2.GetX() + d1.GetY()*d2.GetY() + d1.GetZ()*d2.GetZ());
}

void
test_rotation(const I3Direction &d1, const I3Direction &d2)
{
	I3Direction rotated = RotateToZenith(d1, d2);
	// ENSURE_EQUAL(d1.CalcTheta()-M_PI, -d1.GetZenith());
	// ENSURE_EQUAL(M_PI-d1.CalcPhi(), -d1.GetAzimuth());
	
	ENSURE_EQUAL(Angle(d1, d2)/I3Units::degree, Angle(I3Direction(0,0,-1), rotated)/I3Units::degree);
}

TEST(InTrackSystemIsReversed)
{
	test_rotation(I3Direction(0, 0, -1), I3Direction(0, 0, -1));
	test_rotation(I3Direction(0, -1, 0), I3Direction(0, -1, 0.1));
	// I3Particle track;
	// track.SetPos(0, 0, 0);
	// track.SetDir(0, -1, 0);
	// track.SetShape(I3Particle::InfiniteTrack);
	// 
	// I3Direction odir(0, -1, 0.1);
	// printf("%f %f %f (%f deg)\n", odir.GetX(), odir.GetY(), odir.GetZ(), Angle(track.GetDir(), odir)/I3Units::degree);
	// // I3Direction rotated = I3Calculator::InTrackSystem(track, odir);
	// I3Direction rotated = RotateToZenith(track.GetDir(), odir);
	// printf("%f %f %f (%f deg)\n", rotated.GetX(), rotated.GetY(), rotated.GetZ(), Angle(I3Direction(0,0,-1), rotated)/I3Units::degree);
	// // std::cout << rotated.GetZ() << std::endl;
}
