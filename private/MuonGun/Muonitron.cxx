
#include <MuonGun/Muonitron.h>
#include <MuonGun/CompactTrack.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/I3Calculator.h>
#include <simclasses/I3MMCTrack.h>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

I3_MODULE(Muonitron);

Muonitron::Muonitron(const I3Context &ctx) : I3Module(ctx), spropagator_(new I3MuonGun::MuonPropagator("ice"))
{
	AddParameter("Depths", "Propagate muons to these vertical depths (in meters water-equivalent)", depths_);
	AddParameter("MMC", "Instance of I3PropagatorServiceBase", propagator_);
	AddParameter("CylinderHeight", "Height of the MMC cylinder", 0.);
	
	AddOutBox("OutBox");
}

void
Muonitron::Configure()
{
	GetParameter("Depths", depths_);
	GetParameter("MMC", propagator_);
	GetParameter("CylinderHeight", cyl_length_);
	
	if (depths_.size() == 0)
		log_fatal("You must specify at least one vertical depth!");
	if (!propagator_)
		log_fatal("No MMC propagator configured!");
}

// Convert to a coordinate system where the zenith is given by the given direction
// rather than I3Direction(0,0,-1)
I3Direction
Muonitron::RotateToZenith(const I3Direction &direction, const I3Direction &dir)
{
	I3Direction p(dir);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Position
Muonitron::RotateToZenith(const I3Direction &direction, const I3Position &pos)
{
	I3Position p(pos);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Position
Muonitron::Impact(const I3Particle &p)
{
	// Vector pointing from the origin to the anchor point,
	// projected along the track
	double l = (p.GetPos().GetX()*p.GetDir().GetX()
	          + p.GetPos().GetY()*p.GetDir().GetY()
	          + p.GetPos().GetZ()*p.GetDir().GetZ());
	// Shift anchor point down the track to the closest approach
	// to the origin
	return I3Position(p.GetPos().GetX() - l*p.GetDir().GetX(),
	                  p.GetPos().GetY() - l*p.GetDir().GetY(),
	                  p.GetPos().GetZ() - l*p.GetDir().GetZ());
}

I3Particle
Muonitron::RotateToZenith(const I3Particle &reference, const I3Particle &part)
{
	I3Particle p(part);
	p.SetDir(RotateToZenith(reference.GetDir(), p.GetDir()));
	// Force reference axis to pass through the origin
	I3Position impact = Impact(reference);
	I3Position anchor = I3Position(p.GetPos().GetX()-impact.GetX(),
	                               p.GetPos().GetY()-impact.GetY(),
	                               p.GetPos().GetZ()-impact.GetZ());
	p.SetPos(RotateToZenith(reference.GetDir(), anchor));
	return p;
}

// Find the distance to the surface from a point at depth d
double
Muonitron::GetOverburden(double zenith, double d, double r)
{
	double ct = cos(zenith);
	return sqrt(2*r*d + ct*ct*(r-d)*(r-d) - d*d) - (r-d)*ct;
}

// Transform a detector-centered zenith angle to an earth-centered zenith angle
double
Muonitron::GetGeocentricZenith(double zenith, double d, double r)
{
	double p = GetOverburden(zenith, d, r);
	return atan2(p*sin(zenith), p*cos(zenith) + (r-d));
}


double
Muonitron::GetSurfaceZenith(double zenith, double d, double r)
{
	return zenith - GetGeocentricZenith(zenith, d, r);
}

// Propagate the given muon a set distance, and set its energy and x-y position at
// the observation surface. If the muon decayed before reaching the surface, return false.
bool
Muonitron::PropagateTrack(I3Particle &target, double slant_depth)
{
	std::vector<I3Particle> secondaries;
	
	// Move the new track up in z so that it enters the simulation volume at the desired slant depth 
	double dz = slant_depth + cyl_length_/2.;
	I3Position pos = target.GetPos();
	pos.SetZ(dz);
	target.SetPos(pos);
	I3MMCTrackPtr mmctrack = propagator_->Propagate(target, secondaries);
	if (mmctrack && mmctrack->GetEi() > 0) {
		target.SetPos(mmctrack->GetXi(), mmctrack->GetYi(), 0);
		target.SetTime(mmctrack->GetTi());
		target.SetEnergy(mmctrack->GetEi());
		return true;
	} else {
		// Remove muons that didn't make it
		return false;
	}
}

bool
Muonitron::PropagateTrackSimple(I3Particle &target, double slant_depth)
{	
	target = spropagator_->propagate(target, slant_depth);
	return target.GetEnergy() > 0;
}

void
Muonitron::DAQ(I3FramePtr frame)
{
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	if (!mctree)
		log_fatal("No MCTree!");
	
	I3MCTree::iterator it = mctree->begin();
	const I3Particle &primary = *it;
	
	std::list<I3Particle> tracks;
	for (it++; it != mctree->end(); it++)
		if (it->GetType() == I3Particle::MuPlus || it->GetType() == I3Particle::MuMinus)
			tracks.push_back(RotateToZenith(primary, *it));
	
	// std::cout << "Surface: " << (tracks.size()) << " muons" << std::endl;
	
	I3MuonGun::TrackBundlePtr bundle = boost::make_shared<I3MuonGun::TrackBundle>();
	double traveled = 0;
	BOOST_FOREACH(double vdepth, depths_) {
		// Convert each vertical depth (in meters water-equivalent) to a slant depth in ice,
		// subtracting the portion the muon has already traveled through
		double dx = GetOverburden(primary.GetDir().GetZenith(), vdepth/IceDensity) - traveled;
		std::vector<I3MuonGun::CompactTrack> deep_tracks;
		for (std::list<I3Particle>::iterator pit = tracks.begin(); pit != tracks.end(); ) {
			// if (PropagateTrack(*pit, dx)) {
			if (PropagateTrackSimple(*pit, dx)) {	
				deep_tracks.push_back(I3MuonGun::CompactTrack(*pit));
				pit++;
			} else {
				pit = tracks.erase(pit);
			}
		}
		
		traveled += dx;
		
		// std::cout << "zenith: " << primary.GetDir().GetZenith() << " ice depth: " << vdepth/IceDensity << std::endl;
		// std::cout << vdepth << " mwe (" << traveled << " slant): " << (deep_tracks.size()) << " muons" << std::endl;
		
		if (deep_tracks.size() > 0)
			(*bundle)[vdepth].swap(deep_tracks); 
	}
	
	frame->Put("MCPrimary", boost::make_shared<I3Particle>(primary));
	frame->Put("Tracks", bundle);
	PushFrame(frame);
}
