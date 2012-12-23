/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * @version $Revision$
 * @date $Date$
 */

#ifndef I3MUONGUN_TRACKBINNER_H_INCLUDED
#define I3MUONGUN_TRACKBINNER_H_INCLUDED

#include <MuonGun/histogram.h>
#include <MuonGun/CompactTrack.h>

namespace I3MuonGun {

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

};

#endif