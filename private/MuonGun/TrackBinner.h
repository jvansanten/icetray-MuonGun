/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_TRACKBINNER_H_INCLUDED
#define I3MUONGUN_TRACKBINNER_H_INCLUDED

#include <MuonGun/histogram.h>
#include <MuonGun/CompactTrack.h>
#include <dataclasses/physics/I3MCTree.h>

namespace I3MuonGun {

typedef I3Map<I3ParticleID, double> I3ParticleIDMap;
I3_POINTER_TYPEDEFS(I3ParticleIDMap);

/**
 * @brief A utility class for filling muon bundles into histograms
 */
class TrackBinner {
public:
	TrackBinner(double mindepth, double maxdepth, unsigned steps);
	void Consume(boost::shared_ptr<const TrackBundle> tracks, double energy, double zenith, double weight);

	boost::shared_ptr<histogram::histogram<2> > primary_;
	boost::shared_ptr<histogram::histogram<3> > multiplicity_;
	boost::shared_ptr<histogram::histogram<4> > radius_;
	boost::shared_ptr<histogram::histogram<5> > energy_;
};

class NeutrinoBinner {

public:
	NeutrinoBinner();
	void Consume(boost::shared_ptr<const TrackBundle> tracks,
	    I3MCTreeConstPtr tree, I3ParticleIDMapConstPtr particle_weights, double weight);
	    // I3MCTreeConstPtr tree, double weight);
	typedef histogram::histogram<3> hist;
	typedef boost::shared_ptr<hist> histptr;
	typedef std::map<int, std::vector<histptr> > histmap;
	histmap histograms_;
	
};

};

#endif
