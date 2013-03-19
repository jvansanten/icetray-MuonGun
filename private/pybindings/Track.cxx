
#include "MuonGun/Track.h"
#include "simclasses/I3MMCTrack.h"
#include "icetray/python/list_indexing_suite.hpp"

void
register_Track()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<Track, TrackPtr, bases<I3Particle> >("Track")
	    .def("get_energy", (double (Track::*)(double) const)&Track::GetEnergy)
	    .def("harvest", &Track::Harvest)
	    .staticmethod("harvest")
	;
	
	class_<std::list<Track>, boost::shared_ptr<std::list<Track> > >("TrackList")
	    .def(list_indexing_suite<std::list<Track> >())
	;
}
