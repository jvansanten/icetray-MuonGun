
import numpy
import copy
import sys
from optparse import Option, OptionValueError

# hobo versions of IceTray enumerations
class I3Units:
	GeV = 1.
	TeV = 1e3
	PeV = 1e6
	EeV = 1e9

class ParticleType:
	PPlus       =   14
	He4Nucleus  =  402
	N14Nucleus  = 1407
	Al27Nucleus = 2713
	Fe56Nucleus = 5626

class GenerationProbability(object):
	"""
	A probability distribution from which MC events are drawn.
	"""
	def __init__(self, emin, emax, nevents=1, particle_type=None):
		self.emin = emin
		self.emax = emax
		self.nevents = nevents
		self.particle_type = particle_type
	
	def __imul__(self, factor):
		self.nevents *= factor
		return self
	
	def __mul__(self, factor):
		t = copy.copy(self)
		t *= factor
		return t
	
	def __rmul__(self, factor):
		return self*factor
	
	def __idiv__(self, factor):
		self.nevents /= factor
		return self
		
	def __div__(self, factor):
		t = copy.copy(self)
		t /= factor
		return t
		
	def __iadd__(self, other):
		if isinstance(other, type(self)):
			if self.is_compatible(other):
				self.nevents += other.nevents
				return self
			else:
				raise ValueError("This is a different generation spectrum!")
		else:
			raise TypeError("Can't add a %s to this %s" % (type(other).__name__, type(self).__name__))
	
	def __add__(self, other):
		if isinstance(other, type(self)):
			if self.is_compatible(other):
				t = copy.copy(self)
				t += other
				return t
		else:
			return GenerationProbabilityCollection([self, other])
	
	def __call__(self, E, particle_type=None):
		E = numpy.asarray(E)
		if numpy.ndim(E) == 0:
			if (E>=self.emin and E<self.emax and (particle_type is None or particle_type == self.particle_type)):
				return self.nevents*self.generation_probability(E)
			else:
				return 0
		else:
			mask = (E>=self.emin)&(E<self.emax)
			if particle_type is not None:
				mask &= (particle_type==self.particle_type)
			return numpy.where(mask, self.nevents*self.generation_probability(E), 0)
	
	def is_compatible(self, other):
		raise NotImplementedError("should be overridden")
	
	def generation_probability(self, E):
		raise NotImplementedError("should be overridden")
	
	@classmethod
	def option(cls, to_extend=Option):
		"""
		Extend optparse.Option (or a subclass thereof) to understand GenerationProbability
		"""
		
		def check_option(option, opt, value):
			# only expose subclasses of GenerationProbability
			globalvars = dict([(subcls.__name__, subcls) for subcls in cls.__subclasses__()])
			try:
				genprob = eval(value, dict(), globalvars)
			except Exception:
				e = sys.exc_info()[1]
				raise OptionValueError(str(e))
			if isinstance(genprob, GenerationProbability) or isinstance(genprob, CombinedSample):
				return genprob
			else:
				raise OptionValueError("Can't parse a GenerationProbability out of '%s'" % value)
		
		class GenerationProbabilityOption(to_extend):
			TYPES = to_extend.TYPES + (cls.__name__,)
			TYPE_CHECKER = copy.copy(to_extend.TYPE_CHECKER)
			TYPE_CHECKER[cls.__name__] = check_option
		
		return GenerationProbabilityOption
		

class GenerationProbabilityCollection(object):
	"""
	A collection of generation spectra, possibly for different particle types.
	"""
	def __init__(self, spectra):
		from collections import defaultdict
		self.spectra = defaultdict(list)
		for dist in spectra:
			self.spectra[dist.particle_type].append(dist)
		
	def __call__(self, E, particle_type=None):
		if particle_type is None:
			return sum([prob(E) for spectra in self.spectra.values() for prob in spectra])
		else:
			if numpy.ndim(particle_type) == 0:
				return sum([prob(E) for prob in self.spectra[particle_type]])
			else:
				E = numpy.asarray(E)
				count = numpy.zeros(E.shape)
				for ptype in numpy.unique(particle_type):
					mask = particle_type==ptype
					if numpy.any(mask):
						Em = E[mask]
						count[mask] += sum([prob(Em) for prob in self.spectra[ptype]])
			return count
	
	def __imul__(self, factor):
		for spectra in self.spectra.values():
			for prob in spectra:
				prob *= factor
		return self
	
	def __idiv__(self, factor):
		self *= (1./factor)
		return self
	
	def __mul__(self, factor):
		clone = copy.deepcopy(self)
		clone *= factor
		return clone
	
	def __div__(self, factor):
		return self*(1./factor)
	
	def __rmul__(self, factor):
		return self*factor
	
	def __iadd__(self, other):
		if isinstance(other, type(self)):
			for pt, ospectra in other.spectra.items():
				for ospec in ospectra:
					for spec in self.spectra[pt]:
						if spec.is_compatible(ospec):
							spec += ospec
							break
					else:
						self.spectra[pt].append(ospec)
			return self
		else:
			for spec in self.spectra[other.particle_type]:
				if spec.is_compatible(other):
					spec += other
					break
			else:
				self.spectra[other.particle_type].append(other)
			return self
	
	def __add__(self, other):
		t = copy.copy(self)
		t += other
		return t
	
class PowerLaw(GenerationProbability):
	"""
	Power-law spectra are easy.
	"""
	def __init__(self, eslope, *args, **kwargs):
		super(PowerLaw, self).__init__(*args, **kwargs)
		self.eslope = eslope
		self.gen_norm = self.norm(self.emin, self.emax, self.eslope)
	
	@staticmethod
	def norm(emin, emax, eslope):
		if eslope < -1:
			g = eslope+1
			return (emax**g - emin**g)/g
		else:
			return numpy.log(emax/emin)
	
	def __repr__(self):
		return "PowerLaw(%.2f, emin=%.2e, emax=%.2e, nevents=%.2e)" % (self.eslope, self.emin, self.emax, self.nevents)
	
	def generation_probability(self, E):
		return E**(self.eslope)/self.gen_norm
	
	def invert(self, p):
		"""
		Return CDF^{-1}(p)
		"""
		return (p*(self.emax**(self.eslope+1) - self.emin**(self.eslope+1)) + self.emin**(self.eslope+1))**(1./(self.eslope+1))
	
	def is_compatible(self, other):
		if isinstance(other, type(self)):
			return self.emin == other.emin and self.emax == other.emax and self.eslope == other.eslope and self.particle_type == other.particle_type
		else:
			return False

def FiveComponent(nevents, emin, emax, normalization=[10., 5., 3., 2., 1.], gamma=[-2.]*5, spric=True):
	"""
	Special case: 5-component dCORSIKA
	
	:param norm: relative normalizations of the 5 components
	:param gamma: power law index for each component
	:param spric: make lower energy proportional to primary mass
	"""
	pt = [getattr(ParticleType, p) for p in ('PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus')]
	if spric:
		masses = [[1, p/100][p>100] for p in pt]
	else:
		masses = [1.]*5
	# DCORSIKA does this integral in TeV, so we have to do so as well.
	fluxsums = numpy.array([n*PowerLaw.norm(emin*m/I3Units.TeV, emax/I3Units.TeV, g) for m, g, n in zip(masses, gamma, normalization)])
	nshower = nevents*fluxsums/fluxsums.sum()
	return GenerationProbabilityCollection([PowerLaw(g, emin*m, emax, nevents=n, particle_type=p) for m, g, n, p in zip(masses, gamma, nshower, pt)])

class Hoerandel(GenerationProbability):
	"""
	Special case: DCORSIKA RANPRI=2 mode generates a "natural" spectrum
	"""
	# Ripped from dCORSIKA source (with RANPRI=2)
	z = numpy.arange(26)+1
	# Round atomic masses like dCORSIKA does
	mass_number = numpy.round([1.00797, 4.0026, 6.939, 9.0122, 10.811, 12.0112, 14.0067, 15.9994, 18.9984, 20.183, 22.9898, 24.312, 26.9815, 28.086, 30.984, 32.064, 35.453, 39.948, 39.102, 40.08, 44.956, 47.9, 50.942, 51.996, 54.938, 55.847])
	flux = numpy.array([0.0873, 0.0571, 0.00208, 0.000474, 0.000895, 0.0106, 0.00235, 0.0157, 0.000328, 0.0046, 0.000754, 0.00801, 0.00115, 0.00796, 0.00027, 0.00229, 0.000294, 0.000836, 0.000536, 0.00147, 0.000304, 0.00113, 0.000631, 0.00136, 0.00135, 0.0204])
	gamma = numpy.array([2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64, 2.65, 2.7, 2.64, 2.61, 2.63, 2.67, 2.46, 2.59])
	assert(z.size == mass_number.size)
	def __init__(self, z, dslope=0, *args, **kwargs):
		super(Hoerandel, self).__init__(emin, emax, nevents)
		if z < 1 or z > 26:
			raise ValueError("Unknown nuclear charge %d" % z)
		self.z = z
		zi = int(z)-1
		self.gamma = self.gamma[zi] + dslope
		self.gen_norm = self.fluxsum(self.emin, self.emax, z, self.gamma)/self.primary_probabilities(self.emin, dslope)[zi]
	
	def generation_probability(self, E):
		return self.fluxdiff(E, self.z, self.gamma)/self.gen_norm
	
	@classmethod
	def primary_probabilities(cls, emin, dslope=0):
		# This calculation is done in TeV (it matters!)
		emin = emin*cls.mass_number/I3Units.TeV
		# Generation spectrum includes the slope change
		gamma = cls.gamma + dslope
		inf_fluxsum = cls.flux*(emin**(1-gamma))/(gamma-1)
		return inf_fluxsum/inf_fluxsum.sum()
	
	@staticmethod	
	def fluxdiff(e, z, gamma, delta_gamma=2.1, eps_cutoff=1.9, E_knee=4.49*I3Units.PeV):
		"""
		Differential (unnormalized) Hoerandel flux
		"""
		return e**(-gamma)*(1+(e/(E_knee*z))**eps_cutoff)**(-delta_gamma/eps_cutoff)
	
	@staticmethod
	def fluxsum(emin, emax, z, gamma, delta_gamma=2.1, eps_cutoff=1.9, E_knee=4.49*I3Units.PeV):
		"""
		Integral Hoerandel flux
		"""
		# the Gauss hypergeometric function. whee!
		from scipy.special import hyp2f1
		antideriv = lambda e: ((e**(1-gamma))/(1-gamma))*hyp2f1(delta_gamma/eps_cutoff, (1-gamma)/eps_cutoff, (1-gamma)/eps_cutoff+1, -(e/(E_knee*z))**eps_cutoff)
		return antideriv(emax) - antideriv(emin)

class EnergyWeight(object):
	def __init__(self, target_flux, generation_spectrum):
		self.target_flux = target_flux
		self.generation_spectrum = generation_spectrum
	def __call__(self, E, zenith):
		return self.target_flux(E)/self.generation_spectrum(E, zenith)

