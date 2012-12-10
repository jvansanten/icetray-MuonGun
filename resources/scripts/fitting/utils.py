
import tables, numpy, dashi

def histload(hdf, where):
	"""
	dashi.histload with optional normalization
	"""
	h = dashi.histload(hdf, where)
	node = hdf.getNode(where)
	if 'count' in node._v_attrs:
		h /= node._v_attrs['count']
	print 'norm: %.1f' % node._v_attrs['count']
	return h

def points(self, differential=False):
	sp = dashi.scatterpoints.points2d()
	sp.x = self.bincenters
	sp.xerr = self.binwidths/2.
	sp.y = self.bincontent.copy()
	sp.yerr = self.binerror
	if differential:
		sp.y /= self.binwidths
		sp.yerr /= self.binwidths
	return sp
dashi.histogram.hist1d.points = points

def load_group(fname, group='energy'):
	with tables.openFile(fname) as hdf:
		import operator
		elements = 'He4Nucleus', 'Fe56Nucleus', 'N14Nucleus', 'Al27Nucleus', 'PPlus'
		h = histload(hdf, '/%s/%s' % (elements[0], group))
		for e in elements[1:]:
		    h += histload(hdf, '/%s/%s' % (e, group))
	return h

def load_espec(fname, single=True, bias=50):
	he = load_group(fname, 'energy')

	if single:
		# select single muons and sum over all radii
		he = he[:,:,1,:,:].project([0,1,3])
		eaxis = 2
	else:
		eaxis = 4
		
	# normalize
	norm = he._h_bincontent.sum(axis=eaxis).reshape(he._h_bincontent.shape[:-1] + (1,))
	he._h_bincontent /= norm
	he._h_squaredweights /= norm*norm
	# convert to differential
	shape = (1,)*(he.bincontent.ndim-1) + (he.bincontent.shape[-1],)
	norm = numpy.diff(he._h_binedges[-1][1:-1]).reshape(shape)
	vrange = he._h_visiblerange
	he._h_bincontent[vrange] /= norm
	he._h_squaredweights[vrange] /= norm*norm
	
	# convert energies to log-log
	he._h_bincontent[:] = numpy.log(he._h_bincontent) + bias
	he._h_binedges[-1][1:-1] = numpy.log(he._h_binedges[-1][1:-1])
	
	# convert zenith angles to cos(zenith), reversing the angular axis in the process
	he._h_binedges[0][1:-1] = numpy.cos(he._h_binedges[0][1:-1][::-1])
	# reverse through a temporary to avoid overwriting bits we want to read later
	rev = he._h_bincontent[::-1,:,:].copy()
	he._h_bincontent[:] = rev
	rev = he._h_squaredweights[::-1,:,:].copy()
	he._h_squaredweights[:] = rev
	
	# zero out non-finite weights
	mask = numpy.logical_not(numpy.isfinite(he._h_bincontent) & numpy.isfinite(he._h_squaredweights))
	he._h_bincontent[mask] = 0
	he._h_squaredweights[mask] = 0
	
	return he

def load_radial_distribution(fname, bias=50):
	h = load_group(fname, 'radius')
	#normalize
	norm = h._h_bincontent.sum(axis=3).reshape(h._h_bincontent.shape[:-1] + (1,))
	h._h_bincontent /= norm
	h._h_squaredweights /= norm*norm

	# convert to differential
	shape = (1,)*(h.bincontent.ndim-1) + (h.bincontent.shape[-1],)
	norm = numpy.diff(h._h_binedges[-1][1:-1]**2).reshape(shape)
	vrange = h._h_visiblerange
	h._h_bincontent[vrange] /= norm
	h._h_squaredweights[vrange] /= norm*norm
	
	# parameterize log(dP/dR^2) as a function of R
	h._h_bincontent[:] = numpy.log(h._h_bincontent) + bias
	
	# convert zenith angles to cos(zenith), reversing the angular axis in the process
	h._h_binedges[0][1:-1] = numpy.cos(h._h_binedges[0][1:-1][::-1])
	# reverse through a temporary to avoid overwriting bits we want to read later
	rev = h._h_bincontent[::-1,:,:].copy()
	h._h_bincontent[:] = rev
	rev = h._h_squaredweights[::-1,:,:].copy()
	h._h_squaredweights[:] = rev
	
	# zero out non-finite weights
	mask = numpy.logical_not(numpy.isfinite(h._h_bincontent) & numpy.isfinite(h._h_squaredweights))
	h._h_bincontent[mask] = 0
	h._h_squaredweights[mask] = 0
	
	return h

def pad_knots(knots, order=2):
	"""
	Pad knots out for full support at the boundaries
	"""
	pre = knots[0] - (knots[1]-knots[0])*numpy.arange(order, 0, -1)
	post = knots[-1] + (knots[-1]-knots[-2])*numpy.arange(1, order+1)
	return numpy.concatenate((pre, knots, post))

def colorize(sequence, cmap='jet'):
	from matplotlib.colors import Normalize
	from matplotlib import cm
	norm = Normalize(vmin=0, vmax=len(sequence)-1)
	cmap = getattr(cm, cmap)
	for i, s in enumerate(sequence):
		yield cmap(norm(i)), s

def plot_energy_slice(h, spline, slice_=(3,10,1,1), **kwargs):
	import pylab
	import dashi; dashi.visual()
	from icecube.photospline.glam.glam import grideval
	
	idx = slice_ + (slice(None),)
	sub = h[idx]
	idx = tuple([s-1 for s in slice_]) + (slice(None),)
	
	coords = []
	for i, s in enumerate(idx):
		axis = h._h_bincenters[i][s]
		try:
			len(axis)
			#axis = numpy.linspace(axis[0], axis[-1], 101)
			axis = numpy.linspace(0, 20, 101)
			# print axis
		except TypeError:
			axis = [axis]
		coords.append(axis)
		deg = (180*numpy.arccos(coords[0][0])/numpy.pi)
	sub.scatter(**kwargs)
	pylab.plot(coords[-1], grideval(spline, coords).flatten(), **kwargs)

def plot_radial_slice(h, spline, slice_=(3,10,1), **kwargs):
	import pylab
	import dashi; dashi.visual()
	from icecube.photospline.glam.glam import grideval
	
	idx = slice_ + (slice(None),)
	sub = h[idx]
	idx = tuple([s-1 for s in slice_]) + (slice(None),)
	
	coords = []
	for i, s in enumerate(idx):
		axis = h._h_bincenters[i][s]
		try:
			len(axis)
			#axis = numpy.linspace(axis[0], axis[-1], 101)
			axis = numpy.linspace(0, 250, 1001)
			# print axis
		except TypeError:
			axis = [axis]
		coords.append(axis)
	sub.scatter(**kwargs)
	print grideval(spline, coords).flatten()
	pylab.plot(coords[-1], grideval(spline, coords).flatten(), **kwargs)
