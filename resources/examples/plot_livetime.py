#!/usr/bin/env python

"""
Use muon flux weights to calculate an effective livetime for combined
CORSIKA samples as a function of energy.
"""

import pylab, numpy
from icecube import MuonGun

surface = MuonGun.Cylinder(1600, 800)
area = numpy.pi**2*surface.radius*(surface.radius+surface.length)

# 1 file of E^-2.6 5-component 3e4-1e9 GeV (3:2:1:1:1)
soft = 4e5*MuonGun.corsika_genprob('CascadeOptimized5Comp')
# 1 file of E^-2 5-component 6e2-1e11 GeV (10:5:3:2:1)
hard = 2.5e6*MuonGun.corsika_genprob('Standard5Comp')
# In order to compare to "unweighted" CORSIKA, turn the Hoerandel flux
# into a probability (since we happen to know the integral)
areanorm = 0.131475115*area
# 1 file of natural-spectrum ("unweighted") CORSIKA
unweighted = (2.5e8/areanorm)*MuonGun.corsika_genprob('Hoerandel5')

def get_weight(weighter, energy, zenith=numpy.pi/8, z=0):
	shape = energy.shape
	x = surface.radius*numpy.ones(shape)
	y = numpy.zeros(shape)
	z = z*numpy.ones(shape)
	azimuth = numpy.zeros(shape)
	zenith = zenith*numpy.ones(shape)
	multiplicity = numpy.ones(shape, dtype=numpy.uint32)
	mmax = multiplicity.max()
	e = numpy.zeros(shape + (mmax,), dtype=numpy.float32)
	e[:,0] = energy
	r = numpy.zeros(shape + (mmax,), dtype=numpy.float32)
	return weighter(x, y, z, zenith, azimuth, multiplicity, e, r)

e = numpy.logspace(1, 7, 101)
pylab.figure()
target = MuonGun.load_model('Hoerandel5_atmod12_SIBYLL')
generators = [
	('200k $E^{-2}$ 5-component CORSIKA', 2e5*hard),
	('100k $E^{-2.6}$ 5-component CORSIKA', 1e5*soft),
	('200k unweighted CORSIKA', 2e5*unweighted),
	('total', 2e5*hard + 1e5*soft + 2e5*unweighted),
]

annum = 365*24*3600
for label, generator in generators:
	weighter = MuonGun.WeightCalculator(surface, target, generator)
	pylab.plot(e, 1./(get_weight(weighter, e)*annum), label=label)

pylab.loglog()
pylab.legend(loc='best', prop=dict(size='small'))
pylab.ylabel('Single-muon livetime [years]')
pylab.xlabel('Muon energy at sampling surface [GeV]')
pylab.grid()

pylab.show()
