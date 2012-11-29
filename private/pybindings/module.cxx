
#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>

namespace bp = boost::python;

#include <MuonGun/CompactTrack.h>
#include <MuonGun/Muonitron.h>
#include <MuonGun/I3MuonGun.h>

#include <icetray/python/dataclass_suite.hpp>

void
register_Muonitron()
{
	bp::def("overburden", &Muonitron::GetOverburden, (bp::arg("z"), bp::arg("d")=OriginDepth, bp::arg("r")=SurfaceRadius));
	bp::def("impact", &Muonitron::Impact);
	
	{
		using namespace I3MuonGun;
		bp::class_<Surface, boost::noncopyable>("Surface", bp::no_init)
			.def("intersection", &Surface::GetIntersection)
		;
		
		bp::class_<Cylinder, bp::bases<Surface> >("Cylinder", bp::init<double, double>())
		;
		
		bp::class_<Sphere, bp::bases<Surface> >("Sphere", bp::init<double, double>())
		;
	}

	bp::def("rotate_to_zenith", (I3Particle (*)(const I3Particle&, const I3Particle&))&Muonitron::RotateToZenith);
}

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

void register_histogram();
void register_TrackBinner();
void register_MuonPropagator();
void register_RadialDistribution();
void register_CanCan();

I3_PYTHON_MODULE(MuonGun)
{
	load_project("libMuonGun", false);
	
	register_histogram();
	register_TrackBinner();
	register_MuonPropagator();
	register_extras();
	register_Muonitron();
	register_CompactTrack();
	register_RadialDistribution();
	register_CanCan();
	
}

