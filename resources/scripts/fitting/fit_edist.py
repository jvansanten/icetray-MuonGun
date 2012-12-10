#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--single", dest="single", action="store_true", default=False)

opts, args = parser.parse_args()

import numpy, tables, dashi
from icecube.photospline import spglam as glam

from utils import load_espec, pad_knots

fname = "/data/uwa/jvansanten/projects/2012/muongun/corsika/SIBYLL/Hoerandel5/atmod_12.hdf5"

bias = 50
hes = load_espec(fname, opts.single, bias)

extents = [(e[1], e[-2]) for e in hes._h_binedges]

weights = 1./(numpy.log(numpy.exp(hes.bincontent - bias) + numpy.sqrt(hes._h_squaredweights[hes._h_visiblerange])) - (hes.bincontent - bias))
weights[numpy.logical_not(numpy.isfinite(weights))] = 0
# weights[0,:,:] = 0 # ingore the horizon bin

knots = [
	pad_knots(numpy.linspace(0, 1, 11), 2),  # cos(theta)
	pad_knots(numpy.linspace(1, 3, 11), 2),  # vertical depth [km]
	pad_knots(numpy.linspace(0, 20, 21), 2), # log(energy) [log(E/GeV)]
]
smooth = weights[weights > 0].mean() # make the smoothing term proportional to the weights
order = [2,2,2]
penalties = {2:[smooth/1e3, smooth/1e4, smooth/1e2]}    # Penalize curvature 

if not opts.single:
	knots = knots[0:2] + [
		pad_knots(numpy.linspace(0, numpy.sqrt(100), 3)**2), # bundle multiplicity
		pad_knots(numpy.linspace(0, numpy.sqrt(100), 3)**2), # distance to shower axis
		
	] + [knots[-1]]
	order = order[0:2] + [2,2] + [order[-1]]
	penalties[2] = penalties[2][0:2] + [smooth/1e4, smooth/1e4] + [penalties[2][-1]]

spline = glam.fit(hes.bincontent,weights,hes._h_bincenters,knots,order,smooth,penalties=penalties)
spline.bias = bias

from icecube.photospline import splinefitstable

from icecube.photospline.glam.glam import grideval
import pylab
from utils import colorize
dashi.visual()

for di in xrange(3, 19, 3):
	mi = 1
	ri = 1
	pylab.figure()
	for c, zi in colorize(range(1, 11), 'jet'):
		if opts.single:
			sub = hes[zi,di,:]
			idx = (zi-1, di-1, slice(None))
		else:
			sub = hes[zi,di,mi,ri,:]
			idx = (zi-1, di-1, mi-1, ri-1, slice(None))
		coords = []
		for i, s in enumerate(idx):
			axis = hes._h_bincenters[i][s]
			try:
				len(axis)
				#axis = numpy.linspace(axis[0], axis[-1], 101)
				axis = numpy.linspace(0, 20, 101)
				# print axis
			except TypeError:
				axis = [axis]
			coords.append(axis)
		deg = (180*numpy.arccos(coords[0][0])/numpy.pi)
		sub.scatter(color=c, label='%d deg' % deg)
		print deg, sub.bincontent.sum()
		# print coords
		pylab.plot(coords[-1], grideval(spline, coords).flatten(), color=c)
	pylab.legend()
	title = '%d m' % (coords[1][0]*1000)
	if not opts.single:
		title += ', N~=%.1f, r~=%.1f m' % (hes._h_bincenters[2][mi-1], hes._h_bincenters[3][ri-1])
	pylab.title(title)

pylab.show()