/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_TRACK_H_INCLUDED
#define I3MUONGUN_TRACK_H_INCLUDED

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>
#include <boost/tuple/tuple.hpp>

I3_FORWARD_DECLARATION(I3MMCTrack);
typedef I3Vector<I3MMCTrack> I3MMCTrackList;

namespace I3MuonGun {

/**
 * A particle that moves in a straight line, and loses energy as it does so.
 */
class Track : public I3Particle {
public:
	Track() {};
	Track(const I3MMCTrack &mmctrack,
	    const I3MCTree::sibling_iterator &secondaries_begin,
	    const I3MCTree::sibling_iterator &secondaries_end);
	
	/**
	 * Get the energy of the particle at the given down-track distance,
	 * assuming that the "continuous" contribution to the energy loss
	 * is constant between measurement points.
	 *
	 * @param[in] length distance from the track origin
	 * @returns an energy if 0 <= length < range; otherwise 0
	 */
	double GetEnergy(double length) const;
	I3Position GetPos(double length) const;
	double GetTime(double length) const;
	
	// Un-hide overridden base class methods
	using I3Particle::GetEnergy;
	using I3Particle::GetPos;
	using I3Particle::GetTime;
	
	/**
	 * Extract energy losses from frame objects.
	 *
	 * Find the stochastic energy losses associated with each I3MMCTrack
	 * in the I3MCTree, and store them together in a Track.
	 */
	static std::list<Track> Harvest(const I3MCTree &, const I3MMCTrackList &);

private:
	/**
	 * A point at which the absolute energy of the particle is known
	 */
	struct Checkpoint {
		double length, energy;
		unsigned offset;
		Checkpoint(double l, double e=0., unsigned o=0)
		    : length(l), energy(e), offset(o) {}
	};
	/**
	 * The sum of stochastic energy losses since the last checkpoint
	 */
	struct LossSum {
		double length, energy;
		LossSum(double l, double e=0.)
		    : length(l), energy(e) {}
	};
	std::vector<Checkpoint> checkpoints_;
	std::vector<LossSum> losses_;
};

}

#endif
