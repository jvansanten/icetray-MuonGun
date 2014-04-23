
# a collection of cosmic-ray fluxes

from weighting import I3Units, ParticleType
import numpy
import operator

class CompiledFlux(object):
	"""
	An efficient pre-compiled form of a multi-component flux. For single-element evalutions
	this is ~2 times faster than switching on the primary type with an if statement; for 1e5
	samples it is 2000 times faster than operating on masked slices for each primary type.
	"""
	def __init__(self, expr):
		import numexpr
		self.expr = numexpr.NumExpr(expr)
	def __call__(self, E, ptype):
		return self.expr(E, ptype)
	
	@classmethod
	def build_lookup(cls, mapping, var='ptype', default=0.):
		"""
		Build an expression equivalent to a lookup table
		"""
		if len(mapping) > 0:
			return 'where(%s==%s, %s, %s)' % (var, mapping[0][0], mapping[0][1], cls.build_lookup(mapping[1:], var, default))
		else:
			return str(default)

class Hoerandel(CompiledFlux):
	"""
	All-particle spectrum (up to iron) after Hoerandel_, as implemented
	in dCORSIKA.
	
	.. _Hoerandel: http://dx.doi.org/10.1016/S0927-6505(02)00198-6
	"""
	delta_gamma = 2.1
	eps_cutoff = 1.9
	E_knee = 4.49*I3Units.PeV
	z = "where(ptype > 100, ptype%100, 1)"
	knee = "(1+(E/(%(E_knee)s*%(z)s))**%(eps_cutoff)s)**(-%(delta_gamma)s/%(eps_cutoff)s)" % locals()
	ptypes = [14, 402, 703,  904, 1105, 1206, 1407, 1608, 1909, 2010, 2311, 2412, 2713, 2814, 3115, 3216, 3517, 4018, 3919, 4020, 4521, 4822, 5123, 5224, 5525, 5626]
	def __init__(self):
		gamma = numpy.array([2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64, 2.65, 2.7, 2.64, 2.61, 2.63, 2.67, 2.46, 2.59])
		flux = numpy.array([0.0873, 0.0571, 0.00208, 0.000474, 0.000895, 0.0106, 0.00235, 0.0157, 0.000328, 0.0046, 0.000754, 0.00801, 0.00115, 0.00796, 0.00027, 0.00229, 0.000294, 0.000836, 0.000536, 0.00147, 0.000304, 0.00113, 0.000631, 0.00136, 0.00135, 0.0204])
		flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
		expr = "%(flux)s*E**(-%(gamma)s)*%(knee)s" % dict(gamma=self.build_lookup(list(zip(self.ptypes, gamma))), flux=self.build_lookup(list(zip(self.ptypes, flux))), knee=self.knee)
		CompiledFlux.__init__(self, expr)

class Hoerandel5(Hoerandel):
	"""
	Hoerandel_ with only 5 components, after Becherini_ et al. (also the same as Arne_ Schoenwald's version)
	
	.. _Hoerandel: http://dx.doi.org/10.1016/S0927-6505(02)00198-6
	.. _Arne: http://www.ifh.de/~arnes/Forever/Hoerandel_Plots/
	.. _Becherini: http://dx.doi.org/10.1016/j.astropartphys.2005.10.005
	"""
	ptypes = [getattr(ParticleType, p) for p in ('PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus')]
	def __init__(self):
		gamma = numpy.array([2.71, 2.64, 2.68, 2.67, 2.58])
		flux = numpy.array([8.73e-2, 5.71e-2, 3.24e-2, 3.16e-2, 2.18e-2])
		flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
		expr = "%(flux)s*E**(-%(gamma)s)*%(knee)s" % dict(gamma=self.build_lookup(list(zip(self.ptypes, gamma))), flux=self.build_lookup(list(zip(self.ptypes, flux))), knee=self.knee)
		CompiledFlux.__init__(self, expr)

class Glasstetter(CompiledFlux):
	"""
	Broken power laws for proton and iron components, with spectral
	indices from Glasstetter_ and normalization from Hoerandel_ (via Eike Middell)
	
	FIXME: this predicts a flux that is a factor of ~2 larger than Hoerandel
	
	.. _Glasstetter: http://adsabs.harvard.edu/abs/1999ICRC....1..222G
	.. _Hoerandel: http://dx.doi.org/10.1016/S0927-6505(02)00198-6
	"""
	ptypes = [getattr(ParticleType, p) for p in ('PPlus', 'Fe56Nucleus')]
	def __init__(self):
		gamma1 = self.build_lookup(list(zip(self.ptypes, [2.67, 2.69])))
		gamma2 = self.build_lookup(list(zip(self.ptypes, [3.39, 3.10])))
		flux = self.build_lookup(list(zip(self.ptypes, [17e3, 9e3])))
		knee = self.build_lookup(list(zip(self.ptypes, [4.1*I3Units.PeV, 1.1e2*I3Units.PeV])))
		CompiledFlux.__init__(self, "%(flux)s*where(E < %(knee)s, E**-%(gamma1)s, %(knee)s**(%(gamma2)s-%(gamma1)s)*E**-%(gamma2)s)" % locals())

class GaisserHillas(CompiledFlux):
	"""
	Spectral fits from an `internal report`_ (also on the arXiv) by Tom Gaisser_.
	
	.. _`internal report`: icecube/201102004-v2
	.. _Gaisser: http://arxiv.org/abs/1111.6675v2
	"""
	ptypes = [getattr(ParticleType, p) for p in ('PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus')]
	def get_expression(self, flux, gamma, rigidity):
		z = "where(ptype > 100, ptype%100, 1)"
		return "%(flux)s*E**(-%(gamma)s)*exp(-E/(%(rigidity)s*%(z)s))" % locals()
	def get_flux(self):
		return [[7860., 3550., 2200., 1430., 2120.]]
	def get_gamma(self):
		return [[2.66, 2.58, 2.63, 2.67, 2.63]]
	def get_rigidity(self):
		return [4*I3Units.PeV]
	def __init__(self):
		flux = [self.build_lookup(list(zip(self.ptypes, f))) for f in self.get_flux()]
		gamma = [self.build_lookup(list(zip(self.ptypes, g))) for g in self.get_gamma()]
		rigidity = self.get_rigidity()
		CompiledFlux.__init__(self, "+".join([self.get_expression(f, g, r) for f, g, r in zip(flux, gamma, rigidity)]))

class GaisserH3a(GaisserHillas):
	"""
	The model H3a with a mixed extra-galactic population (Fig. 2)
	has all iron at the highest energy and would represent a
	scenario in which the cutoff is not an effect of energy loss
	in propagation over cosmic distances through the CMB but is
	instead just the extra-galactic accelerators reaching their
	highest energy.
	"""
	def get_flux(self):
		return super(GaisserH3a, self).get_flux() + [[20]*2 + [13.4]*3, [1.7]*2 + [1.14]*3]
	def get_gamma(self):
		return super(GaisserH3a, self).get_gamma() + [[2.4]*5, [2.4]*5]
	def get_rigidity(self):
		return super(GaisserH3a, self).get_rigidity() + [30*I3Units.PeV, 2*I3Units.EeV]

class GaisserH4a(GaisserH3a):
	"""
	In the model H4a, on the other hand, the extra-galactic component
	is assumed to be all protons.
	"""
	def get_flux(self):
		return super(GaisserH4a, self).get_flux()[:-1] + [[200]]
	def get_gamma(self):
		return super(GaisserH4a, self).get_gamma()[:-1] + [[2.6]]
	def get_rigidity(self):
		return super(GaisserH4a, self).get_rigidity()[:-1] + [60*I3Units.EeV]
