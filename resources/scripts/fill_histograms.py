#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray

import sys
infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', filenamelist=infiles)

import numpy
from icecube.icetray import I3Units
class HoerandelWeight(object):
	gamma = numpy.array([2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64, 2.65, 2.7, 2.64, 2.61, 2.63, 2.67, 2.46, 2.59])
	flux = numpy.array([0.0873, 0.0571, 0.00208, 0.000474, 0.000895, 0.0106, 0.00235, 0.0157, 0.000328, 0.0046, 0.000754, 0.00801, 0.00115, 0.00796, 0.00027, 0.00229, 0.000294, 0.000836, 0.000536, 0.00147, 0.000304, 0.00113, 0.000631, 0.00136, 0.00135, 0.0204])
	flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
	def __init__(self, emin, emax, eslope, z=1, nevents=1e6*1):
		self.nevents = nevents
		self.emin = emin
		self.emax = emax
		self.eslope = eslope
		self.z = z
		if self.eslope < -1:
			g = self.eslope+1
			self.gen_norm = (self.emax**g - self.emin**g)/g
		else:
			self.gen_norm = numpy.log(self.emax/self.emin)
		self.norm = self.flux[int(z)-1]/self.nevents#/self.fluxsum(self.emin, self.emax, z, self.gamma[int(z)-1])
		
	def generation_probability(self, E):
		return E**(self.eslope)/self.gen_norm
	
	def fluxdiff(self, e, z, gamma, delta_gamma=2.1, eps_cutoff=1.9, E_knee=4.49*I3Units.PeV):
		"""
		Differential (unnormalized) Hoerandel flux
		"""
		return e**(-gamma)*(1+(e/(E_knee*z))**eps_cutoff)**(-delta_gamma/eps_cutoff)
	
	def fluxsum(self, emin, emax, z, gamma, delta_gamma=2.1, eps_cutoff=1.9, E_knee=4.49*I3Units.PeV):
		"""
		Integral Hoerandel flux
		"""
		# the Gauss hypergeometric function. whee!
		from scipy.special import hyp2f1
		antideriv = lambda e: ((e**(1-gamma))/(1-gamma))*hyp2f1(delta_gamma/eps_cutoff, (1-gamma)/eps_cutoff, (1-gamma)/eps_cutoff+1, -(e/(E_knee*z))**eps_cutoff)
		return antideriv(emax) - antideriv(emin)
		
	def __call__(self, E):
		#return 1/(self.nevents*self.generation_probability(E))
		return self.norm*self.fluxdiff(E, self.z, self.gamma[int(self.z)-1])/self.generation_probability(E)

class HoerandelWeight5(HoerandelWeight):
	"""
	Hoerandel with only 5 components, after Becherini et al. (also the same as Arne Schoenwald's version)
	"""
	gamma = numpy.nan*numpy.zeros(26)
	flux  = numpy.nan*numpy.zeros(26)
	gamma[0]  = 2.71 # H
	gamma[1]  = 2.64 # He
	gamma[6]  = 2.58 # N
	gamma[12] = 2.67 # Al
	gamma[25] = 2.58 # Fe
	flux[0]  = 8.73e-2 # H
	flux[1]  = 5.71e-2 # He
	flux[6]  = 3.24e-2 # N
	flux[12] = 3.16e-2 # Al
	flux[25] = 2.18e-2 # Fe
	flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
	def __init__(self, *args, **kwargs):
		HoerandelWeight.__init__(self, *args, **kwargs)
		if numpy.isnan(self.norm):
			raise ValueError("I can't account for nuclei with charge %d" % self.z)

import dashi, numpy, tables
from icecube import MuonGun
class Filla(icetray.I3Module):
	def __init__(self, ctx):
		icetray.I3Module.__init__(self, ctx)
		self.AddOutBox("OutBox")
	
	def Configure(self):
		from collections import defaultdict
		
		depthbins = numpy.linspace(1.0, 5.0, 9)
		depthbins -= numpy.diff(depthbins)[0]/2.
		zenbins = numpy.arccos(numpy.linspace(1, 0, 11))
		zenbins_fine = numpy.arccos(numpy.linspace(1, 0, 101))
		multbins = numpy.array([1, 2, 3, 4, 10, 20, 40, 100], dtype=float)
		rbins = numpy.array([0, 5, 10, 15, 25, 45], dtype=float)
		
		self.primary = dashi.histogram.hist2d((zenbins, numpy.logspace(2, 11, 101)))
		self.multiplicity = dashi.histogram.histogram(3, (zenbins_fine, depthbins, numpy.arange(1, 100)))
		self.radius = dashi.histogram.histogram(4, (zenbins, depthbins, multbins, numpy.linspace(0, numpy.sqrt(250), 101)**2))
		self.energy = dashi.histogram.histogram(5, (zenbins, depthbins, multbins, rbins, numpy.logspace(0, 6, 101)))
		
		
		self.multiplicity_slices = tuple([tuple([self.multiplicity[i+1,j+1,:] for j in xrange(len(depthbins))]) for i in xrange(len(zenbins_fine))])
		self.radius_slices = tuple([tuple([self.radius[i+1,j+1,:,:] for j in xrange(len(depthbins))]) for i in xrange(len(zenbins))])
		self.energy_slices = tuple([tuple([self.energy[i+1,j+1,:,:,:] for j in xrange(len(depthbins))]) for i in xrange(len(zenbins))])
		
		
		self.depthbins = depthbins
		self.zenbins = zenbins
		self.zenbins_fine = zenbins_fine
		
		self.weighter = None
		
		import os
		if os.path.exists(outfile):
			os.unlink(outfile)
			
		self.nevents = 0
		
	def DAQ(self, frame):
		primary = frame['MCPrimary']
		#if primary.type != primary.PPlus:
		#	return
		# print primary
		if self.weighter is None:
			if 'CorsikaWeightDict' in frame:
				# generated by I3CORSIKAReader
				wm = frame['CorsikaWeightDict']
				if primary.type == primary.PPlus:
					z = 1
				else:
					z = int(primary.type)%100
				self.weighter = HoerandelWeight5(wm['EnergyPrimaryMin'], wm['EnergyPrimaryMax'], wm['PrimarySpectralIndex'], z, wm['NEvents'])
			elif 'CorsikaWeightMap' in frame:
				# generated by I3CORSIKAWeightModule
				wm = frame['CorsikaWeightMap']
				if wm['Weight'] != 1.0:
					raise ValueError("Can't deal with weighted DCorsika")
				timescale = wm['TimeScale']
				area = wm['AreaSum']
				self.weighter = lambda energy: 1/(timescale*area)
		weight = self.weighter(primary.energy)
		#print weight
		
		#frame['Tracks'].items()
		#return
		
		zenith = numpy.asarray([primary.dir.zenith])
		zi = max(numpy.searchsorted(self.zenbins, zenith)[0] - 1, 0)
		zif = max(numpy.searchsorted(self.zenbins_fine, zenith)[0] - 1, 0)
		
		self.primary.fill((numpy.asarray([zenith]), numpy.asarray([primary.energy])), weights=numpy.asarray([weight]))
		
		multiplicity=self.multiplicity_slices[zif]
		radius=self.radius_slices[zi]
		energy=self.energy_slices[zi]
		
		for di, (depth, tracks) in enumerate(frame['Tracks'].iteritems()):
			#print zi, zenith[0], di, depth
			kmwe = depth/I3Units.km
			mult = len(tracks)
			values = numpy.asarray([(mult, p.radius, p.energy) for p in tracks])
			
			one = numpy.ones(1)
			#self.multiplicity[di+1,zi+1,:].fill(mult*one, weights=weight*one)
			multiplicity[di].fill(mult*one, weights=weight*one)
			
			w = (weight/len(tracks))*numpy.ones(len(tracks))
			#self.radius[di+1, zi+1, :, :].fill(values[:,:2], weights=w)
			#self.energy[di+1, zi+1, :, :, :].fill(values, weights=w)
			radius[di].fill(values[:,:2], weights=w)
			energy[di].fill(values, weights=w)
		
		self.nevents += 1
		if self.nevents % 100 == 0:
			print '%d events' % self.nevents
		
		self.PushFrame(frame)
		
	def Finish(self):
		for i in xrange(len(self.zenbins_fine)):
			for j in xrange(len(self.depthbins)):
				self.multiplicity._h_bincontent[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_bincontent
				self.multiplicity._h_squaredweights[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_squaredweights
		for i in xrange(len(self.zenbins)):
			for j in xrange(len(self.depthbins)):
				self.radius._h_bincontent[i+1,j+1,:,:] += self.radius_slices[i][j]._h_bincontent
				self.radius._h_squaredweights[i+1,j+1,:,:] += self.radius_slices[i][j]._h_squaredweights
				self.energy._h_bincontent[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_bincontent
				self.energy._h_squaredweights[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_squaredweights
				
		with tables.openFile(outfile, 'w') as hdf:
			dashi.histsave(self.primary, hdf, '/', 'primary')
			dashi.histsave(self.multiplicity, hdf, '/', 'multiplicity')
			dashi.histsave(self.radius, hdf, '/', 'radius')
			dashi.histsave(self.energy, hdf, '/', 'energy')
		
tray.AddModule(Filla, 'filla')

# tray.AddModule('I3Writer', 'writer',
#     Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
#     DropOrphanStreams=[icetray.I3Frame.DAQ],
#     filename=outfile)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
