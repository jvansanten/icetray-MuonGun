
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
	h._h_squaredweights / norm*norm
	
	fig = pylab.figure()
	fig.subplots_adjust(left=0.15)
	di = 5
	for color, zi in colorize(range(1, 11)):
		sub = h[zi,di,:]
		zenith = h._h_bincenters[0][zi-1]
		ct = numpy.cos(zenith)
		depth = h._h_bincenters[1][di-1]
		sub.scatter(differential=False, color=color, label='%d deg' % (180*zenith/numpy.pi))
		if normed:
			norm = 1
		else:
			norm = flux(depth, ct)
		pylab.plot(sub.bincenters, sub.binwidths*norm*numpy.array([efit(depth, ct, 1, 0, e) for e in sub.bincenters]), color=color)
	if log:
		pylab.loglog()
	else:
		pylab.semilogx()
	if normed:
		if log:
			pylab.ylim((1e-15, 1e-1))
		pylab.ylabel('$dP/dE\,[1/GeV]$')
	else:
		if log:
			pylab.ylim((1e-25, 1e-2))
		pylab.ylabel(log_fluxlabel(sub.binedges))
	pylab.xlabel('Muon energy [GeV]')
	pylab.legend(prop=dict(size='medium'), loc='best', ncol=2)
	pylab.grid()
	pylab.title('Single-muon energy distribution at %d m, various zenii' % (depth*1000))
	
	fig = pylab.figure()
	fig.subplots_adjust(left=0.15)
	zi = 3
	for color, di in colorize(range(5, 19, 2)):
		sub = h[zi,di,:]
		zenith = h._h_bincenters[0][zi-1]
		ct = numpy.cos(zenith)
		depth = h._h_bincenters[1][di-1]
		sub.scatter(differential=False, color=color, label='%d m' % (numpy.round(depth*1000)))
		if normed:
			norm = 1
		else:
			norm = flux(depth, ct)
		pylab.plot(sub.bincenters, sub.binwidths*norm*numpy.array([efit(depth, ct, 1, 0, e) for e in sub.bincenters]), color=color)
	if log:
		pylab.loglog()
	else:
		pylab.semilogx()
	if normed:
		if log:
			pylab.ylim((1e-10, 1e-1))
	else:
		if log:
			pylab.ylim((1e-12, 1e-4))
		pylab.ylabel(log_fluxlabel(sub.binedges))
		
	pylab.xlabel('Muon energy [GeV]')
	pylab.legend(prop=dict(size='medium'), loc='best', ncol=2)
	pylab.grid()
	pylab.title('Single-muon energy distribution at %d deg zenith, various depths' % (180*zenith/numpy.pi))

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
		pylab.plot(sub.bincenters, numpy.array([flux(depth, ct) for ct in coszen]), color=color)
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
		curve = numpy.array([flux(depth, ct) for ct in coszen])
		chi2 = ((sp.y-curve)/sp.yerr)**2
		chi2 = chi2[numpy.isfinite(chi2)].sum()
		print chi2
		sp.y /= curve
		sp.yerr /= curve
		mask = numpy.isnan(sp.y)
		sp.y[mask] = 1
		sp.yerr[mask] = 1
		sp.scatter(color=color, label='%d m' % (numpy.round(depth*1000)))
	pylab.ylim((0.5, 1.5))
	pylab.ylabel('ratio MC/fit')
	pylab.grid()
	pylab.xlabel('Zenith angle [degree]')
	
	return ax1, ax2
	
	