
#ifndef I3MUONGUN_TRACKBINNER_H_INCLUDED
#define I3MUONGUN_TRACKBINNER_H_INCLUDED

#include <MuonGun/histogram.h>
#include <MuonGun/CompactTrack.h>

namespace I3MuonGun {
	
// typedef I3Map<double, std::vector<CompactTrack> > CompactTrackSeriesMap;
// I3_POINTER_TYPEDEFS(CompactTrackSeriesMap);

class TrackBinner {
public:
	TrackBinner(double mindepth, double maxdepth, unsigned steps);
	void Consume(boost::shared_ptr<const TrackBundle> tracks, double energy, double zenith, double weight);

	boost::shared_ptr<histogram<2> > primary_;
	boost::shared_ptr<histogram<3> > multiplicity_;
	boost::shared_ptr<histogram<4> > radius_;
	boost::shared_ptr<histogram<5> > energy_;
};

};

#endif