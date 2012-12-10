#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--single", dest="single", action="store_true", default=False)

opts, args = parser.parse_args()

import numpy, tables, dashi
from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable

from utils import load_radial_distribution, pad_knots

fname = "/data/uwa/jvansanten/projects/2012/muongun/corsika/SIBYLL/Hoerandel5/atmod_12.hdf5"

bias = 50
h = load_radial_distribution(fname, bias)

extents = [(e[1], e[-2]) for e in h._h_binedges]

weights = 1./(numpy.log(numpy.exp(h.bincontent - bias) + numpy.sqrt(h._h_squaredweights[h._h_visiblerange])) - (h.bincontent - bias))
weights[numpy.logical_not(numpy.isfinite(weights))] = 0
weights[:,:,0,:] = 0 # ignore the single-track bin

knots = [
	pad_knots(numpy.linspace(0, 1, 11), 2),  # cos(theta)
	pad_knots(numpy.linspace(1, 3, 11), 2),  # vertical depth [km]
	pad_knots(numpy.linspace(0, numpy.sqrt(100), 6)**2, 2), # bundle multiplicity
	pad_knots(numpy.linspace(0, numpy.sqrt(250), 21)**2, 2), # radius^2 [m^2]
]
smooth = weights[weights > 0].mean()/1e6 # make the smoothing term proportional to the weights
order = [2,2,2,2]
penalties = {2:[smooth/1e3, smooth/1e4, smooth/1e4, smooth/1e2]}    # Penalize curvature 

spline = glam.fit(h.bincontent,weights,h._h_bincenters,knots,order,smooth,penalties=penalties)
spline.bias = bias

import os
os.unlink('Hoerandel5_atmod12_SIBYLL.radius.fits')
splinefitstable.write(spline, 'Hoerandel5_atmod12_SIBYLL.radius.fits')
