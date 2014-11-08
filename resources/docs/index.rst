.. 
.. copyright  (C) 2013
.. Jakob van Santen <vansanten@wisc.edu>
.. and The Icecube Collaboration http://www.icecube.wisc.edu
.. 
.. $Id$
.. 
.. @version $Revision$
.. @date $LastChangedDate$
.. @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy$

.. highlight:: python

.. _MuonGun-main:

MuonGun
=======

MuonGun is a toolkit for efficiently simulating the flux of atmospheric muons deep under the ice.

Rationale
^^^^^^^^^

In diffuse neutrino searches it is almost always necessary to estimate the
background due to atmospheric muons from simulation. In IceCube this has
typically been done directly, by simulating air showers to ground level with
CORSIKA, propagating the muons in the shower through the firn and ice with
MMC to a cylindrical sampling surface surrounding the detector, and then
weighting the simulated events to an assumed cosmic-ray flux. Though this
method offers the highest possible precision available from the chosen
simulation software, it suffers from two key inefficiencies. First, since
the simulation starts with cosmic-ray primaries rather than in-ice muons,
one has only loose control over the characteristics of the muon bundles that
actually reach the detector. For example, if one were interested only in
single muons with a few TeV of energy, one would spend quite a lot of time
simulating both showers whose muons never reached the detector those that
result in high-multiplicity bundles. Second, the direct approach makes it
necessary to repeat the entire simulation chain in order to change aspects
of the air shower simulation such as atmospheric profile or hadronic model.

An alternative approach is to de-couple the air shower simulation and muon
propagation from the remainder of the simulation by constructing a
parametrization of the muon flux under the ice and drawing muon bundles from
the parameterized distribution. This allows one to generate specific bundle
configurations and weight them properly, and also to re-weight existing
simulated events to a muon flux associated with different assumptions about
interactions in the atmosphere.

The parametric approach is used heavily by ANTARES in the form of
their MUPAGE_ event generator. The work described here is an attempt to
apply the technique, described in a paper by `Becherini et al`_, to IceCube
simulation.

Flux parametrization
^^^^^^^^^^^^^^^^^^^^

MuonGun works by drawing samples from a parametrization of the atmospheric muon
flux as a function of vertical depth, zenith angle, multiplicity, energy, and
for bundles, the distance of each muon from the shower axis. These variables
(the same ones used in a parametrization by `Becherini et al`_) completely
describe a muon bundle under the following assumptions:

1. The flux of cosmic-ray primaries that reach the atmosphere and their 
   daughter muons is independent of the azimuthal arrival direction.
2. The bundle is azimuthally symmetric around the shower axis.
3. All muons in the bundle are perfectly parallel to each other and to the 
   shower axis.
4. The shower front has no curvature.

The first assumption is violated by deflection in the Earth's magnetic field,
but the same approximation is used whenever dCORSIKA showers are randomized in
azimuth before being fed in to IceTray (i.e. nearly always). The remaining
three approximations are important only if it is possible to measure very
detailed properties of the bundle structure over relatively short (~1 km)
observation baselines.

The parameterizations are made by propagating muons from dCORSIKA simulation
filling their properties into a set of histograms
(resources/scripts/propagate_and_fill.py). The distributions in these
histograms are then fitted with tensor-product B-spline surfaces using
photospline (resources/fitting/fit.sh).

Generating muons
^^^^^^^^^^^^^^^^

Muon bundle generators are instances of the :cpp:class:`Generator` class. There are currently 3 generators implemented:

.. cpp:class:: StaticSurfaceInjector
	
	Injects muon bundles on a fixed sampling surface with depth and zenith 
	angle distributions given by a flux model and an energy distribution 
	given by a broken power law. This is more or less identical to MUPAGE_.

.. cpp:class:: EnergyDependentSurfaceInjector
	
	Similar to :cpp:class:`StaticSurfaceInjector`, but scales the sampling 
	surface with energy. This allows for much higher effective livetimes 
	for low energy events in the center of the detector where veto 
	techniques can be most effective.

.. cpp:class:: Floodlight
	
	A demonstration generator that simply illuminates its sampling surface 
	from all directions with single muons.

To generate muon bundles, you must first configure a :cpp:class:`Generator` and
pass it to the provided :py:func:`GenerateBundles` segment. For example, to
generate single muons from the flux of Hoerandel_::
	
	from icecube.icetray import I3Units
	from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
	from icecube.MuonGun.segments import GenerateBundles

	# Use Hoerandel as a template for generating muons
	model = load_model('Hoerandel5_atmod12_SIBYLL')
	# Generate only single muons, no bundles
	model.flux.max_multiplicity = 1
	# Center the sampling surface on the barycenter of IC79 strings
	surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(31.25, 19.64, 0))
	# Draw energies from an E^-2 power law broken at 1 TeV, from 10 TeV to 10 PeV
	spectrum = OffsetPowerLaw(2, 1*I3Units.TeV, 10*I3Units.TeV, 10*I3Units.PeV)
	# Set up the generator. This gets stored in a special frame for later reference
	generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

	tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=10000, GCDFile=gcd)

Weighting
^^^^^^^^^

In order to weight generated muon bundles to a flux, you need to know both the
distribution of events in the flux model (parameterized in the tables provided
with MuonGun) and the distribution of events that you generated (calculated by
the Generator, stored in a special "S" frame at the beginning of every
generated file). Given those, you can calculate weights either within IceTray
using the :cpp:class:`WeightCalculatorModule` I3Module or from a standalone
Python script using the :cpp:class:`WeightCalculator` class.

First, you should collect the generators for all of the files you plan to use::
	
	def harvest_generators(infiles):
		"""
		Harvest serialized generator configurations from a set of I3 files.
		"""
		from icecube.icetray.i3logging import log_info as log
		generator = None
		for fname in infiles:
			f = dataio.I3File(fname)
			fr = f.pop_frame(icetray.I3Frame.Stream('S'))
			f.close()
			if fr is not None:
				for k in fr.keys():
					v = fr[k]
					if isinstance(v, MuonGun.GenerationProbability):
						log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
						if generator is None:
							generator = v
						else:
							generator += v
		return generator

The generators can simply be added together, or multiplied by an integer to
represent a larger number of identically-configured generators. You can pass
this combined generator to :cpp:class:`WeightCalculatorModule` to calculate a
weight appropriate for the combined set of files::
	
	model = MuonGun.load_model(model)
	generator = harvest_generators(infiles)
	tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model,
	    Surface=generator.surface, Generator=generator)

This will put an I3Double called "MuonWeight" into the frame that represents a
weight in events per second. Alternatively, you can use the provided
:ref:`tableio-main` converter to write the parameters needed for the weight
calculation to a table::
	
	from icecube.hdfwriter import I3HDFWriter
	tray.AddSegment(I3HDFWriter, 'scribe',
	    Output=outfile,
	    Keys=[dict(key='I3MCTree', name='BundleParameters',
	             converter=MuonGun.converters.MuonBundleConverter(1, generator.surface))],
	    Types=[],
	    SubEventStreams=['nullsplit'],
	)

and then use the standalone :cpp:class:`WeightCalculator` class to calculate a
weight::
	
	model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
	generator = harvest_generators(infiles)
	weighter = MuonGun.WeightCalculator(generator.surface, model, generator)
	
	with tables.openFile(outfile) as hdf:
		axis = hdf.root.MCPrimary.read()
		bundle = hdf.root.BundleParameters.read()
		weights = weighter(axis['x'], axis['y'], axis['z'], axis['zenith'], axis['azimuth'],
		    bundle['multiplicity'], bundle['energy'], bundle['radius'])

.. note:: The weighter will only be able to accept Numpy arrays if you have `boost::numpy`_ installed. If you do not have `boost::numpy`_ it will simply be exposed as a scalar function.

.. _`boost::numpy`: https://github.com/martwo/BoostNumpy/

API
^^^

The class structure is documented in the autogenerated Doxygen_ docs.

.. _Doxygen: ../../doxygen/MuonGun/index.html

.. _`Becherini et al`: http://dx.doi.org/10.1016/j.astropartphys.2005.10.005
.. _MUPAGE: http://dx.doi.org/10.1016/j.cpc.2008.07.014
.. _Hoerandel: http://dx.doi.org/10.1016/S0927-6505(02)00198-6

Utilities: recovering the energies of muons at depth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Muons are ranged particles, so both their position and energy are a function of
time. The position is easy enough to represent in an :cpp:class:`I3Particle`;
by dint of special relativity, the position be approximated by speed-of-light
displacement along the muon's initial direction from an arbitrary reference
point. Information about a muon's true energy at any given moment, however, is
scattered in multiple locations in each frame. These are:

- The initial energy at the reference point, stored in the :cpp:class:`I3Particle`
  that represents the muon itself
- The energy the muon has when it enters and exits the MMC simulation volume
  (typically a 1600 x 800 meter upright cylinder), stored in the :cpp:class:`I3MMCTrack`
  corresponding to the muon.
- The sizes of the stochastic energy losses inside the MMC simulation volume. These
  are stored as :cpp:class:`I3Particle` attached as daughters of the :cpp:class:`I3Particle`
  represending the muon in the :cpp:class:`I3MCTree`

.. py:currentmodule:: icecube.MuonGun

While it is in principle straightforward to recover the energy of a simulated
muon at any point along its path inside the MMC simulation volume from these
scattered data structures, it is also extremely tedious. MuonGun includes
a utility class, :py:class:`Track`, to automate this task.

.. py:class:: Track
	
	A subclass of I3Particle that includes the particle's energy losses
	
	.. classmethod:: harvest(frame)
		
		Assemble a collection of :py:class:`Track` from the
		:cpp:class:`MMCTrackList` and :cpp:class:`I3MCTree` in *frame*
	
	.. method:: get_energy(displacement)
		
		Get the energy of the muon *displacement* meters from its starting position

This class can be used, for example, to find the total energy losses
(stochastic and continuous) of all muons within some volume::
	from icecube import MuonGun, simclasses
	
	# A surface approximating the actual detector (make it smaller if you only care e.g. about DeepCore)
	surface = MuonGun.Cylinder(1000,500)

	edep = 0
	for track in MuonGun.Track.harvest(frame['I3MCTree'], frame['MMCTrackList']):
		# Find distance to entrance and exit from sampling volume	
		intersections = surface.intersection(track.pos, track.dir)
		# Get the corresponding energies
		e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)
		# Accumulate
		edep +=  (e0-e1)

There is also a convenience function, :py:func:`muons_at_surface`, that uses
:py:class:`Track` to "shift" all the muons in a frame forward so that their
positions, times, and energies correspond to their intersections with a given
surface.

.. py:function:: muons_at_surface(frame, surface)
	
	:param frame: an I3Frame
	:param surface: a :py:class:`Surface`
	:returns: a list of I3Particles with positions, times, and energies that
	          correspond to their intersection with *surface*. Muons that range
	          out before reaching *surface* are not included.

This can be used to quickly obtain the true multiplicity of a muon bundle when
it enters the detector::
	surface = MuonGun.Cylinder(1000,500)
	multiplicity = len(MuonGun.muons_at_surface(frame, surface))

Weighting CORSIKA simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In order to make the parameterization it was necessary to weight CORSIKA
air-shower simulation to a number of different flux models. The Python code to
do this is included in MuonGun for reference. We calculate a weight that turns
a number of simulated events into a rate of muon bundles crossing a sampling
surface by dividing the differential flux by the differential number of
simulated events:

.. math::

	w \,\, [\rm s^{-1}] = \frac{d\Phi/dE}{dN/dEdAd\Omega} \left[ \frac{\rm GeV^{-1} s^{-1} m^{-2} sr^{-1}}{ \rm GeV^{-1} m^{-2} sr^{-1}} \right]

.. py:currentmodule:: icecube.MuonGun.fluxes

For the numerator, three parameterizations of the primary cosmic ray flux are provided:

.. autoclass:: Hoerandel5

.. autoclass:: GaisserH3a

.. autoclass:: GaisserH4a

.. py:currentmodule:: icecube.MuonGun.weighting

The denominator is in principle easy to calculate, as it consists of a series
of power laws for each primary nucleus. However, the bookkeeping becomes quite
involved once simulation sets with different settings are combined. The
:py:func:`FiveComponent` factory function simplifies the calculation by
creating a :py:class:`GenerationProbabilityCollection` object that can
calculate a differential number of events for each primary type.

.. autofunction:: FiveComponent

The notion of combining simulation sets can be expressed by addition and multiplication of the appropriate :py:class:`GenerationProbability` objects. For instance, to create a normalization term for two different 5-component sets, one with 100k files and one with 27k::
	
	>>> low = FiveComponent(10000000, 600, 1e5, [5, 2.25, 1.1, 1.2, 1.0], [-2.65, -2.6, -2.6, -2.6, -2.6], )
	<icecube.MuonGun.weighting.GenerationProbabilityCollection object at 0x1099667d0>
	>>> high = FiveComponent(30000, 1e5, 1e11, [3, 2.25, 1.1, 1.2, 1.0], [-2.]*5,  spric=False)
	<icecube.MuonGun.weighting.GenerationProbabilityCollection object at 0x109966990>
	>>> norm = 1e5*low + 2.7e4*high
	<icecube.MuonGun.weighting.GenerationProbabilityCollection at 0x109966410>

Since all of these functions operate on Numpy arrays, they can easily be used
to calculate weights from HDF5 files. As an example, suppose there is an HDF5
file with tables /MCPrimary (for the primary nucleus) and /CorsikaWeightMap
(for the sampling surface area). The normalization calculated above can be used
to make a weight::
	
	energy = hdf.root.MCPrimary.cols.energy[:]
	ptype = hdf.root.MCPrimary.cols.type[:]
	area = hdf.root.CorsikaWeightMap.cols.AreaSum[:]
	flux = MuonGun.Hoerandel5()
	
	weight = area*flux(energy, ptype)/norm(energy, ptype)




