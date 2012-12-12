
import pylab
import numpy
import dashi; dashi.visual()
from utils import load_group, colorize, points

fluxlabel = '$d\Phi/dE\,\,[1/(GeV\,m^2 sr\,s)]$'
def log_fluxlabel(edges):
	d = numpy.diff(numpy.log(edges))[0]
	return '$%.2f \\times d\Phi/d\log10{(E/GeV)}\,\,[1/(m^2 sr\,s)]$' % d

def plot_single_energy(histogram_fname, flux, efit, normed=False, log=True):
	
	h = load_group(histogram_fname, 'energy')[:,:,1,:,:].project([0,1,3])
	if normed:
		# dP/dE
		norm = h._h_bincontent.sum(axis=2).reshape(h._h_bincontent.shape[:-1] + (1,))
	else:
		# dPhi/dE
		norm = numpy.diff(numpy.cos(h._h_binedges[0][::-1])).reshape((h._h_bincontent.shape[0],) + (1,)*(h.ndim-1))	
	
	h._h_bincontent /= norm
	h._h_squaredweights /= norm*norm
	
	fig = pylab.figure()
	fig.subplots_adjust(left=0.15)
	di = 5
	for color, zi in colorize(range(1, 11)):
		sub = h[zi,di,:]
		zenith = h._h_bincenters[0][zi-1]
		ct = numpy.cos(zenith)
		depth = h._h_bincenters[1][di-1]
		if normed:
			d = numpy.diff(numpy.log(sub.binedges))[0]
			norm = 1./d
		else:
			norm = flux(depth, ct, 1)
		(sub/d).scatter(differential=False, color=color, label='%d deg' % (180*zenith/numpy.pi))
		pylab.plot(sub.bincenters, sub.binwidths*norm*numpy.array([efit(depth, ct, 1, 0, e) for e in sub.bincenters]), color=color)
	if log:
		pylab.loglog()
	else:
		pylab.semilogx()
	if normed:
		if log:
			pylab.ylim((1e-12, 1e0))
		pylab.ylabel('$dP/d\log10(E_{\\mu}/GeV) $')
	else:
		if log:
			pylab.ylim((1e-25, 1e-2))
		pylab.ylabel(log_fluxlabel(sub.binedges))
	pylab.xlabel('Muon energy [GeV]')
	pylab.legend(prop=dict(size='medium'), loc='best', ncol=2)
	pylab.grid()
	pylab.title('Single-muon energy distribution at %d m' % (depth*1000))
	
	fig = pylab.figure()
	fig.subplots_adjust(left=0.15)
	zi = 3
	for color, di in colorize(range(5, 19, 2)):
		sub = h[zi,di,:]
		zenith = h._h_bincenters[0][zi-1]
		ct = numpy.cos(zenith)
		depth = h._h_bincenters[1][di-1]
		if normed:
			d = numpy.diff(numpy.log(sub.binedges))[0]
			norm = 1/d
		else:
			norm = flux(depth, ct, 1)
		(sub/d).scatter(differential=False, color=color, label='%d m' % (numpy.round(depth*1000)))
		pylab.plot(sub.bincenters, sub.binwidths*norm*numpy.array([efit(depth, ct, 1, 0, e) for e in sub.bincenters]), color=color)
	if log:
		pylab.loglog()
	else:
		pylab.semilogx()
	if normed:
		if log:
			pylab.ylim((1e-12, 1e0))
		pylab.ylabel('$dP/d\log10(E_{\\mu}/GeV) $')		
	else:
		if log:
			pylab.ylim((1e-12, 1e-4))
		pylab.ylabel(log_fluxlabel(sub.binedges))
		
	pylab.xlabel('Muon energy [GeV]')
	pylab.legend(prop=dict(size='medium'), loc='best', ncol=2)
	pylab.grid()
	pylab.title('Single-muon energy distribution at %d deg zenith' % (180*zenith/numpy.pi))

def plot_single_flux(histogram_fname, flux, log=True):
	h = load_group(histogram_fname, 'multiplicity')[:,:,1]
	norm = numpy.diff(numpy.cos(h._h_binedges[0][::-1])).reshape((h._h_bincontent.shape[0],) + (1,)*(h.ndim-1))	
	h._h_bincontent /= norm
	h._h_squaredweights /= norm*norm
	return plot_single_flux_impl(h, flux, log)
	
def plot_single_flux_impl(h, flux, log=True):
	from matplotlib.gridspec import GridSpec
	from matplotlib.ticker import NullFormatter
	
	h._h_binedges[0] /= (numpy.pi/180)
	
	depths = lambda: range(4, 19, 2)
	
	fig = pylab.figure()
	grid = GridSpec(2,1,height_ratios=[2,1], hspace=0.05)
	ax1 = pylab.subplot(grid[0])
	for color, di in colorize(depths()):
		sub = h[:,di]
		depth = h._h_bincenters[1][di-1]
		sub.scatter(color=color, label='%d m' % (numpy.round(depth*1000)))
		coszen = numpy.cos(numpy.pi*sub.bincenters/180)
		pylab.plot(sub.bincenters, numpy.array([flux(depth, ct, 1) for ct in coszen]), color=color)
	if log:
		pylab.semilogy()
		pylab.ylim((1e-12, 1e-2))
	pylab.ylabel('$\Phi \, [1/(m^2 sr \, s)]$')
	pylab.legend(prop=dict(size='medium'), loc='best', ncol=2)
	pylab.grid()
	pylab.title('Single-muon flux at various depths')
	ax1.get_xaxis().set_major_formatter(NullFormatter())
	
	ax2 = pylab.subplot(grid[1])
	coszen = numpy.cos(numpy.pi*h._h_bincenters[0]/180)
	for color, di in colorize(depths()):
		sp = h[:,di].points()
		depth = h._h_bincenters[1][di-1]
		curve = numpy.array([flux(depth, ct, 1) for ct in coszen])
		chi2 = ((sp.y-curve)/sp.yerr)**2
		chi2 = chi2[numpy.isfinite(chi2)].sum()
		# print chi2
		sp.y /= curve
		sp.yerr /= curve
		mask = numpy.isnan(sp.y)
		sp.y[mask] = 1
		sp.yerr[mask] = 1
		sp.scatter(color=color, label='%d m' % (numpy.round(depth*1000)))
	pylab.ylim((0.95, 1.05))
	pylab.ylabel('ratio MC/fit')
	pylab.grid()
	pylab.xlabel('Zenith angle [degree]')
	
	return ax1, ax2

def plot_bundle_flux(histogram_fname, flux, log=True):
	from utils import load_flux
	return plot_bundle_flux_impl(load_flux(histogram_fname, transform=False), flux, log)
	# h = load_group(histogram_fname, 'multiplicity')
	# norm = numpy.diff(numpy.cos(h._h_binedges[0][::-1])).reshape((h._h_bincontent.shape[0],) + (1,)*(h.ndim-1))	
	# h._h_bincontent /= norm
	# h._h_squaredweights /= norm*norm
	# return plot_bundle_flux_impl(h, flux, log)

def plot_bundle_flux_impl(h, flux, log=True):
	from matplotlib.gridspec import GridSpec
	from matplotlib.ticker import NullFormatter
	
	h._h_binedges[0] /= (numpy.pi/180)
	
	depths = lambda: range(4, 19, 2)
	
	for zi, ylim in zip((2, 45, 60, 75), ((1e-11,1e-2), (1e-12,1e-3), (1e-15,1e-3), (1e-18,1e-4))):
	
		fig = pylab.figure()
		grid = GridSpec(2,1,height_ratios=[2,1], hspace=0.05)
		ax1 = pylab.subplot(grid[0])
		ax2 = pylab.subplot(grid[1])
		
		maxval = -1
		for color, di in colorize(range(4, 19, 2)):
			sub = h[zi,di,:]
			coszen = numpy.cos(numpy.pi*h._h_bincenters[0][zi-1]/180)
			depth = h._h_bincenters[1][di-1]
			pylab.sca(ax1)
			sub.scatter(color=color, label='%d m' % numpy.round(depth*1000))
			curve = numpy.array([flux(depth, coszen, int(m)) for m in sub.bincenters])
			# print sub.bincenters.astype(int)
			pylab.plot(sub.bincenters, curve, color=color)
			if sub.bincontent.max() > maxval:
				maxval = sub.bincontent.max()
			pylab.sca(ax2)
			sp = sub.points()
			sp.y /= curve
			sp.yerr /= curve
			mask = numpy.isnan(sp.y)
			sp.y[mask] = 1
			sp.yerr[mask] = 1
			sp.scatter(color=color)
			
		if log:
			ax1.semilogy()
			maxval = 10**numpy.ceil(numpy.log10(maxval))
			# ax1.set_ylim((maxval/1e16, maxval))
			ax1.set_ylim(ylim)
		ax1.legend(prop=dict(size='medium'), loc='best', ncol=2)
		ax1.set_title('Muon bundle flux at %d deg' % (h._h_bincenters[0][zi-1]))
		ax1.grid()
		ax1.get_xaxis().set_major_formatter(NullFormatter())
		ax1.set_ylabel('$\Phi \, [1/(m^2 sr \, s)]$')
		
		ax2.set_ylim((0.5, 1.5))
		ax2.set_ylabel('ratio MC/fit')
		ax2.grid()
		ax2.set_xlabel('Bundle multiplicity')

def plot_radial_distribution(fname, rdist, log=True):
	from utils import load_radial_distribution
	from matplotlib.gridspec import GridSpec
	from matplotlib.ticker import NullFormatter
	
	h = load_radial_distribution(fname, transform=False)
	h._h_binedges[0] /= (numpy.pi/180)
	
	depths = lambda: range(4, 19, 2)
	
	fig = pylab.figure()
	grid = GridSpec(2,1,height_ratios=[2,1], hspace=0.05)
	ax1 = pylab.subplot(grid[0])
	ax2 = pylab.subplot(grid[1])
	zi = 4
	for color, di in colorize(range(4, 19, 2)):
		m = 2
		sub = h[zi,di,2,:]
		coszen = numpy.cos(numpy.pi*h._h_bincenters[0][zi-1]/180)
		depth = h._h_bincenters[1][di-1]
		
		pylab.sca(ax1)
		sub.scatter(color=color, label='%d m' % numpy.round(depth*1000))
		curve = numpy.array([rdist(depth, coszen, m, r) for r in sub.bincenters])
		pylab.plot(sub.bincenters, curve, color=color)
		
		pylab.sca(ax2)
		sp = sub.points()
		sp.y /= curve
		sp.yerr /= curve
		mask = numpy.isnan(sp.y)
		sp.y[mask] = 1
		sp.yerr[mask] = 1
		sp.scatter(color=color)
	
	ax1.set_xlim((0,50))
	ax2.set_xlim((0,50))
	
	ax1.legend(prop=dict(size='medium'), loc='best', ncol=2)
	ax1.set_title('2-muon radial distribution at %d deg' % (h._h_bincenters[0][zi-1]))
	ax1.grid()
	ax1.get_xaxis().set_major_formatter(NullFormatter())
	ax1.set_ylabel('$dP/dr \, [1/m]$')
		
	ax2.set_ylim((0.95, 1.05))
	ax2.set_ylabel('ratio MC/fit')
	ax2.grid()
	ax2.set_xlabel('Distance to shower axis [m]')
	
	fig = pylab.figure()
	grid = GridSpec(2,1,height_ratios=[2,1], hspace=0.05)
	ax1 = pylab.subplot(grid[0])
	ax2 = pylab.subplot(grid[1])
	di = 10
	for color, zi in colorize(range(1, 8)):
		m = 2
		sub = h[zi,di,2,:]
		coszen = numpy.cos(numpy.pi*h._h_bincenters[0][zi-1]/180)
		depth = h._h_bincenters[1][di-1]
		
		pylab.sca(ax1)
		sub.scatter(color=color, label='%d deg' % (h._h_bincenters[0][zi-1]))
		curve = numpy.array([rdist(depth, coszen, m, r) for r in sub.bincenters])
		pylab.plot(sub.bincenters, curve, color=color)
		
		pylab.sca(ax2)
		sp = sub.points()
		sp.y /= curve
		sp.yerr /= curve
		mask = numpy.isnan(sp.y)
		sp.y[mask] = 1
		sp.yerr[mask] = 1
		sp.scatter(color=color)
	
	ax1.set_xlim((0,50))
	ax2.set_xlim((0,50))
	
	ax1.legend(prop=dict(size='medium'), loc='best', ncol=2)
	ax1.set_title('2-muon radial distribution at %d m' % numpy.round(h._h_bincenters[1][di-1]*1000))
	ax1.grid()
	ax1.get_xaxis().set_major_formatter(NullFormatter())
	ax1.set_ylabel('$dP/dr \, [1/m]$')
		
	ax2.set_ylim((0.95, 1.05))
	ax2.set_ylabel('ratio MC/fit')
	ax2.grid()
	ax2.set_xlabel('Distance to shower axis [m]')
	
	fig = pylab.figure()
	grid = GridSpec(2,1,height_ratios=[2,1], hspace=0.05)
	ax1 = pylab.subplot(grid[0])
	ax2 = pylab.subplot(grid[1])
	di = 10
	zi = 4
	for color, mi in colorize(range(2, 8)):
		sub = h[zi,di,mi,:]
		coszen = numpy.cos(numpy.pi*h._h_bincenters[0][zi-1]/180)
		mult = int(h._h_bincenters[2][mi-1])
		depth = h._h_bincenters[1][di-1]
		
		pylab.sca(ax1)
		sub.scatter(color=color, label='$%d \\leq N_{\\mu} < %d$' % (h._h_binedges[2][mi]+0.5, h._h_binedges[2][mi+1]+0.5))
		curve = numpy.array([rdist(depth, coszen, mult, r) for r in sub.bincenters])
		pylab.plot(sub.bincenters, curve, color=color)
		
		pylab.sca(ax2)
		sp = sub.points()
		sp.y /= curve
		sp.yerr /= curve
		mask = numpy.isnan(sp.y)
		sp.y[mask] = 1
		sp.yerr[mask] = 1
		sp.scatter(color=color)
	
	ax1.set_xlim((0,50))
	ax2.set_xlim((0,50))
	
	ax1.legend(prop=dict(size='medium'), loc='best', ncol=2)
	ax1.set_title('Bundle radial distribution at %d m/%d deg' % (numpy.round(h._h_bincenters[1][di-1]*1000), h._h_bincenters[0][zi-1]))
	ax1.grid()
	ax1.get_xaxis().set_major_formatter(NullFormatter())
	ax1.set_ylabel('$dP/dr \, [1/m]$')
		
	ax2.set_ylim((0.95, 1.05))
	ax2.set_ylabel('ratio MC/fit')
	ax2.grid()
	ax2.set_xlabel('Distance to shower axis [m]')

def plot_bundle_energy_distribution(fname, edist, log=True):
	from utils import load_radial_distribution
	from matplotlib.gridspec import GridSpec
	from matplotlib.ticker import NullFormatter
	
	from utils import load_espec
	h = load_espec(fname, transform=False, single=False)
	h._h_binedges[0] /= (numpy.pi/180)
	
	def plot_edist(label):
		sub = h[zi,di,mi,ri,:]
		radius = h._h_bincenters[3][ri-1]
		coszen = numpy.cos(numpy.pi*h._h_bincenters[0][zi-1]/180)
		mult = int(h._h_bincenters[2][mi-1])
		depth = h._h_bincenters[1][di-1]
		curve = numpy.array([edist(depth, coszen, mult, radius, e) for e in sub.bincenters])
		sub.scatter(color=color, label=label)
		pylab.plot(sub.bincenters, curve, color=color)
		
	def format_axis(title):
		ax1 = pylab.gca()
		ax1.loglog()
		ax1.set_ylim((1e-14, 1e-2))
	
		ax1.legend(prop=dict(size='medium'), loc='best', ncol=1)
		ax1.set_title(title)
		ax1.grid()
		ax1.set_xlabel('Muon energy [GeV] [m]')
		ax1.set_ylabel('$dP/dE_{\\mu} \\, [1/GeV]$')
	
	fig = pylab.figure()
	di, zi, ri, mi = 10, 4, 3, 2
	for color, mi in colorize(range(2, 8)):
		plot_edist('$%d \\leq N_{\\mu} < %d$' % (h._h_binedges[2][mi]+0.5, h._h_binedges[2][mi+1]+0.5))
	format_axis('Bundle energy distribution at %d m/%d deg/%d m off axis' % (numpy.round(h._h_bincenters[1][di-1]*1000), h._h_bincenters[0][zi-1], h._h_bincenters[3][ri-1]))
	
	fig = pylab.figure()
	di, zi, ri, mi = 10, 4, 3, 2
	for color, ri in colorize(range(1, 6)):
		plot_edist('$%d \\, m \\, \\leq r < %.0f \\, m$' % (h._h_binedges[3][ri], h._h_binedges[3][ri+1]))
	format_axis('2-muon energy distribution at %d m/%d deg' % (numpy.round(h._h_bincenters[1][di-1]*1000), h._h_bincenters[0][zi-1]))
	
	fig = pylab.figure()
	di, zi, ri, mi = 10, 4, 3, 2
	for color, zi in colorize(range(1, 8)):
		plot_edist('%d deg' % (h._h_bincenters[0][zi-1]))
	format_axis('2-muon energy distribution at %d m/%d m off axis' % (numpy.round(h._h_bincenters[1][di-1]*1000), h._h_bincenters[3][ri-1]))
	
	fig = pylab.figure()
	di, zi, ri, mi = 10, 4, 3, 2
	for color, di in colorize(range(4, 10)):
		plot_edist('%d m' % numpy.round(h._h_bincenters[1][di-1]*1000))
	format_axis('2-muon energy distribution at %d deg/%d m off axis' % (h._h_bincenters[0][zi-1], h._h_bincenters[3][ri-1]))
	
		