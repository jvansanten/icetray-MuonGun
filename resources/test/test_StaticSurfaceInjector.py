#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, sim_services, simclasses, MuonGun
from I3Tray import I3Tray

from utils import get_propagator
propagator = get_propagator(seed=1337)
if not propagator:
	import sys
	print('Skipping test')
	sys.exit(0)

tray = I3Tray()

tray.AddModule('I3InfiniteSource', 'driver')
tray.AddModule('I3MCEventHeaderGenerator', 'headergen', IncrementEventID=True)
tray.AddService('I3GSLRandomServiceFactory', 'rng', Seed=1337)

surface  = MuonGun.Cylinder(1600, 800)
flux     = MuonGun.BMSSFlux()
flux.min_multiplicity = 1
flux.max_multiplicity = 1
energies = MuonGun.BMSSEnergyDistribution()
radii    = MuonGun.BMSSRadialDistribution()
generator = 10*MuonGun.StaticSurfaceInjector(surface, flux, MuonGun.OffsetPowerLaw(2, 500., 50, 1e6), radii)

tray.AddModule('I3MuonGun::GeneratorModule', 'generator', Generator=generator, Propagator=propagator)
tray.AddModule('I3MuonGun::WeightCalculatorModule', 'weight',
    Flux=flux, EnergyDistribution=energies, RadialDistribution=radii,
    Generator=generator)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
