from icecube.load_pybindings import load_pybindings
import icecube.icetray, icecube.dataclasses # be nice and pull in our dependencies

load_pybindings(__name__,__path__)

def corsika_genprob(atmosphere=12, hadronic_model='sibyll'):
	from os.path import expandvars, join
	basedir = expandvars('$I3_BUILD/MuonGun/resources/scripts/fitting')
	config = 'CascadeOptimized'
	base = join(basedir, '%s_atmod%d_%s' % (config, atmosphere, hadronic_model.upper()))
	return CORSIKAGenerationProbability(Cylinder(1600, 800),
	    SplineFlux(base+'.single_flux.fits', base+'.bundle_flux.fits'),
	    SplineRadialDistribution(base+'.radius.fits'),
	    SplineEnergyDistribution(base+'.single_energy.fits', base+'.bundle_energy.fits'))

