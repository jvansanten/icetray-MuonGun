#!/usr/bin/env python

# specialization for 5-component dCORSIKA

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
from icecube import MuonGun
from I3Tray import I3Tray
from cubicle.weighting import GenerationProbability, PowerLaw, EnergyWeight
from cubicle import fluxes
import numpy
from optparse import OptionParser
parser = OptionParser(option_class=GenerationProbability.option())
parser.add_option("--flux", dest="flux", default="Hoerandel5", choices=("GaisserH3a", "GaisserH4a", "Hoerandel5", "CascadeOptimized"), help="CR flux model to weight to")
parser.add_option("--detcfg", dest="detcfg", type=float, default=None, help="height/diameter ratio of target volume. If unspecified, assume a uniform flux over all solid angle.")
parser.add_option("--mindepth", dest="mindepth", type=float, default=1.0, help="minimum vertical depth to propagate to [%default km]")
parser.add_option("--maxdepth", dest="maxdepth", type=float, default=5.0, help="maximum vertical depth to propagate to [%default km]")
parser.add_option("--steps", dest="steps", type=int, default=9, help="number of vertical depths [%default]")

opts, args = parser.parse_args()
if len(args) < 2:
	parser.error("You must supply at least one input and output file")
infiles, outfile = args[:-1], args[-1]

tray = I3Tray()

from os.path import expandvars
icetray.load('ucr-icetray')
ucr_opts = expandvars('$I3_BUILD/bin/ucr-icetray-ucr ')
for fn in infiles:
	ucr_opts += fn + ' '
ucr_opts += ("-DEPTH=1950 -LENGTH=1600 -RADIUS=800" % locals())
tray.AddModule('I3InfiniteSource', 'driver')
tray.AddModule('I3GeneratorUCR', 'reader', EventsToIssue=int(1e9), UCROpts=ucr_opts)

from utils import MMCFactory, dcorsika_spectra, EnergyWeightCollection, IsotropicWeight, VolumeCorrWeight, Filler, FastFiller, Router
from icecube.MuonGun import MuonPropagator, Crust, Sphere

crust = Crust(MuonPropagator("air", ecut=-1, vcut=5e-2, rho=0.673))
crust.add_layer(Sphere(1948, 6374134), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=0.832))
crust.add_layer(Sphere(1748, 6373934), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005))

tray.AddModule('Muonitron', 'propagator',
    Depths=list(numpy.linspace(opts.mindepth, opts.maxdepth, opts.steps)*I3Units.km),
    Propagator=MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005),
    Crust=crust,
)

ptype = dataclasses.I3Particle.ParticleType
elements = [('H', ptype.PPlus), ('He', ptype.He4Nucleus), ('N', ptype.N14Nucleus), ('Al', ptype.Al27Nucleus), ('Fe', ptype.Fe56Nucleus)]
generated_components = dcorsika_spectra(nevents=2.5e6)

print generated_components

if opts.flux == "CascadeOptimized":
	# Weight to a pseudo-flux that integrates to NEvents/(Area*SolidAngle)
	flux_components = dcorsika_spectra([-2.6]*5, [3., 2., 1., 1., 1.], 3e4, 1e9, 1)
	
	print flux_components
	
	# flux_components = dcorsika_spectra(nevents=1.0)
	r = 800
	l = 1600
	area = numpy.pi**2*r*(r+l)
	# shove the sampling area into the denominator
	generation_spectra = dict([(primary, VolumeCorrWeight(opts.detcfg, f*area)) for (name, primary), f in zip(elements, generated_components)])
	
	fluxes = dict([(primary, f) for (name, primary), f in zip(elements, flux_components)])
	weighter = EnergyWeightCollection(fluxes, generation_spectra)
	tray.AddModule(FastFiller, 'filler',
	    Outfile=outfile, Weight=weighter,
	    MinDepth=opts.mindepth, MaxDepth=opts.maxdepth, DepthSteps=opts.steps)
	tray.AddModule('TrashCan', 'YesWeCan')
	    
else:
	"""
	spectra = dict()
	fluxdict = dict()
	for i, ((name, primary), spectrum) in enumerate(zip(elements, generated_components)):
		if primary == ptype.PPlus:
			z = 1
		else:
			z = int(primary)%100
		flux = getattr(fluxes, opts.flux)(z)
		fluxdict[primary] = flux
		if opts.detcfg:
			spectra[primary] = (VolumeCorrWeight(opts.detcfg, spectrum))
		else:
			spectra[primary] = (IsotropicWeight(spectrum))
	weighter = EnergyWeightCollection(fluxdict, spectra)
	tray.AddModule(FastFiller, 'filler',
	    Outfile=outfile, Weight=weighter,
	    MinDepth=opts.mindepth, MaxDepth=opts.maxdepth, DepthSteps=opts.steps)
	tray.AddModule('TrashCan', 'YesWeCan')	
	
	"""
	
	from collections import defaultdict
	routes = defaultdict(list)
	for i, ((name, primary), spectrum) in enumerate(zip(elements, generated_components)):
		if primary == ptype.PPlus:
			z = 1
		else:
			z = int(primary)%100
		flux = getattr(fluxes, opts.flux)(z)
		if opts.detcfg:
			weighter = EnergyWeight(flux, VolumeCorrWeight(opts.detcfg, spectrum))
		else:
			weighter = EnergyWeight(flux, IsotropicWeight(spectrum))
		name = str(primary)
		tray.AddModule(FastFiller, name,
		    Outfile=outfile+'/'+name, Weight=weighter,
		    MinDepth=opts.mindepth, MaxDepth=opts.maxdepth, DepthSteps=opts.steps)
		routes[primary].append(name)

	tray.AddModule(Router, 'router', Routes=routes)

	tray.AddModule('TrashCan', 'YesWeCan')

	tray.ConnectBoxes('driver', 'OutBox', 'reader')
	tray.ConnectBoxes('reader', 'OutBox', 'propagator')
	tray.ConnectBoxes('propagator', 'OutBox', 'router')
	tray.ConnectBoxes('router', 'OutBox', 'YesWeCan')

	for group in routes.itervalues():
		for name in group:
			tray.ConnectBoxes('router', name, name)
			tray.ConnectBoxes(name, 'OutBox', 'YesWeCan')

tray.Execute()
tray.Finish()
