from icecube.load_pybindings import load_pybindings
import icecube.icetray, icecube.dataclasses # be nice and pull in our dependencies
load_pybindings(__name__,__path__)

from . extruded_polygon import ExtrudedPolygon

def load_model(base):
	from os.path import exists, expandvars, join
	if not exists(base+'.single_flux.fits'):
		basedir=expandvars('$I3_BUILD/MuonGun/resources/tables')
		icecube.icetray.logging.log_info("The table path %s does not exist! Assuming you meant %s..." % (base, join(basedir, base)), unit="MuonGun")
		base = join(basedir, base)
	return BundleModel(SplineFlux(base+'.single_flux.fits', base+'.bundle_flux.fits'),
	    SplineRadialDistribution(base+'.radius.fits'),
	    SplineEnergyDistribution(base+'.single_energy.fits', base+'.bundle_energy.fits'))

from os.path import expandvars
def corsika_genprob(config='CascadeOptimized5Comp', atmosphere=12, hadronic_model='sibyll',
    basedir=expandvars('$I3_BUILD/MuonGun/resources/tables')):
	from os.path import join
	base = join(basedir, '%s_atmod%d_%s' % (config, atmosphere, hadronic_model.upper()))
	return CORSIKAGenerationProbability(Cylinder(1600, 800),
	    SplineFlux(base+'.single_flux.fits', base+'.bundle_flux.fits'),
	    SplineRadialDistribution(base+'.radius.fits'),
	    SplineEnergyDistribution(base+'.single_energy.fits', base+'.bundle_energy.fits'))
del expandvars
