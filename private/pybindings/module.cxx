
#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>

namespace bp = boost::python;

#include <MuonGun/CompactTrack.h>

#include <icetray/python/dataclass_suite.hpp>

void
register_CompactTrack()
{
	using namespace I3MuonGun;
	
	bp::class_<CompactTrack>("CompactTrack")
	    #define PROPS (Energy)(Radius)(Time)(Type)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CompactTrack, PROPS)
	    #undef PROPS
	    .def(bp::dataclass_suite<CompactTrack>())
	;
	
	bp::class_<std::vector<CompactTrack> >("CompactTrackSeries")
	    .def(bp::dataclass_suite<std::vector<CompactTrack> >())
	;
	
	bp::class_<TrackBundle, boost::shared_ptr<TrackBundle>, bp::bases<I3FrameObject> >("CompactTrackSeriesMap")
	    .def(bp::dataclass_suite<TrackBundle>())
	;
	
	register_pointer_conversions<TrackBundle>();
}

void
register_extras()
{
	bp::class_<std::vector<I3Particle::ParticleType> >("I3ParticleTypeSeries")
		.def(bp::dataclass_suite<std::vector<I3Particle::ParticleType> >())
	;
}

void register_RadialDistribution();

I3_PYTHON_MODULE(MuonGun)
{
	load_project("libMuonGun", false);
	
	register_extras();
	register_CompactTrack();
	register_RadialDistribution();
}

