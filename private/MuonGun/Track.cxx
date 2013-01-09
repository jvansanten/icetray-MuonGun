/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/Track.h>
#include <simclasses/I3MMCTrack.h>

#include <boost/foreach.hpp>

namespace I3MuonGun {

Track::Track(const I3MMCTrack &mmctrack,
    const I3MCTree::sibling_iterator &sbegin,
    const I3MCTree::sibling_iterator &send) : I3Particle(mmctrack.GetI3Particle())
{
	// In the beginning, the particle is at its vertex with the given energy
	checkpoints_.push_back(Checkpoint(0., I3Particle::GetEnergy(), 0u));
	if (mmctrack.GetEi() > 0) {
		// Track started outside the MMC volume; we get an extra
		// measurement point (but no stochastics)
		losses_.push_back(LossSum(checkpoints_.back().length, 0.));
		checkpoints_.push_back(Checkpoint((mmctrack.GetTi()-I3Particle::GetTime())/I3Particle::GetSpeed(),
		    mmctrack.GetEi(), losses_.size()));
	}
	
	// Sum energy losses between entry and exit
	double elost = 0;
	BOOST_FOREACH(const I3Particle &p, std::make_pair(sbegin, send)) {
		elost += p.GetEnergy();
		losses_.push_back(LossSum((p.GetTime()-I3Particle::GetTime())/GetSpeed(),
		    elost));
	}
	
	if (mmctrack.GetEf() > 0) {
		// Track made it to the edge of the MMC volume
		losses_.push_back(LossSum(checkpoints_.back().length, elost));
		checkpoints_.push_back(Checkpoint((mmctrack.GetTf()-I3Particle::GetTime())/I3Particle::GetSpeed(),
		    mmctrack.GetEf(), losses_.size()));
		elost = 0.;
	}
	
	losses_.push_back(LossSum(I3Particle::GetLength(), elost));
	checkpoints_.push_back(Checkpoint(I3Particle::GetLength(), 0., losses_.size()));
}

namespace {

template <typename T>
inline bool Sort(const T &a, const T &b)
{
	return a.length < b.length;
}

}

double
Track::GetEnergy(double length) const
{
	if (!std::isfinite(length) || length < 0 || length >= GetLength())
		return 0.;
	// Find an energy checkpoint. The above if() guarantees that
	// that both cp and cp+1 are valid.
	std::vector<Checkpoint>::const_iterator cp =
	    std::lower_bound(checkpoints_.begin(), checkpoints_.end(),
	    Checkpoint(length), Sort<Checkpoint>);
	// Store iterators the mark the records of stochastic losses
	// between the checkpoints
	std::vector<LossSum>::const_iterator l1(losses_.begin()+cp->offset),
	    l2(losses_.begin()+(cp+1)->offset);
	// Find the cumulative energy loss since the last checkpoint
	std::vector<LossSum>::const_iterator ls = std::lower_bound(l1, l2,
	    LossSum(length), Sort<LossSum>);
	
	// Estimate continuous loss rate
	double conti_rate = ((cp+1)->energy - cp->energy - (l2-1)->energy)
	    /((cp+1)->length - cp->length);
	return cp->energy - ls->energy - conti_rate*(length-cp->length);
	
	return 0.;
}

I3Position
Track::GetPos(double length) const
{
	if (!std::isfinite(length) || length < 0 || length >= GetLength())
		return I3Position(NAN, NAN, NAN);
	const I3Position &pos = I3Particle::GetPos();
	const I3Direction &dir = I3Particle::GetDir();
	return I3Position(pos.GetX() + length*dir.GetX(),
	                  pos.GetY() + length*dir.GetY(),
	                  pos.GetZ() + length*dir.GetZ());
}

double
Track::GetTime(double length) const
{
	if (!std::isfinite(length) || length < 0 || length >= GetLength())
		return NAN;
	else
		return I3Particle::GetTime() + length/I3Particle::GetSpeed();
}

namespace {

inline bool
operator!=(const I3MMCTrack &track, const I3Particle &p)
{
	return (track.GetI3Particle().GetMajorID() != p.GetMajorID() ||
	    track.GetI3Particle().GetMinorID() != p.GetMinorID());
}

}

std::list<Track>
Track::Harvest(const I3MCTree &mctree, const I3MMCTrackList &mmctracks)
{
	std::list<Track> tracks;
	I3MCTree::iterator p = mctree.begin();
	BOOST_FOREACH(const I3MMCTrack &mmctrack, mmctracks) {
		// Walk the MCTree for an entry corresponding to this MMCTrack
		while (p != mctree.end() && mmctrack != *p)
			p++;
		if (p != mctree.end()) {
			// Get energy checkpoints from the MMCTrack and stochastic losses
			// from the direct daughters of the corresponding MCTree node
			tracks.push_back(Track(mmctrack, mctree.begin(p), mctree.end(p)));
			// Fast-forward past the secondaries
			p = mctree.end(p);
		}
	}
	
	return tracks;
}

};
