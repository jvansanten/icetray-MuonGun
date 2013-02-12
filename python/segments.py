
from icecube.icetray import traysegment, I3Units
from icecube import dataclasses

import random
from os.path import expandvars
def MakePropagator(radius=800*I3Units.m, length=1600*I3Units.m,
    seed=random.randint(0, (1<<32) - 1), impl='mmc',
    mediadef=expandvars('$I3_BUILD/mmc-icetray/resources/mediadef')):
	"""
	Create a muon propagator service.
	
	:param radius: radius of the target cylinder
	:param length: full height of the target cylinder
	:param impl: if "mmc", use MMC, otherwise use PROPOSAL
	:param mediadef: path to MMC media definition file
	"""
	# Now create the MMC propagators, but first *all* of the options must be set here. 
	# There's no special options added behind the scenes.  This is much more flexible. 
	#  Below are the standard options.  To interpret them see the MMC docs.
	mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
	mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
	mmcOpts += expandvars("-mediadef=%s " % mediadef)
	mmcOpts += "-radius=%d " % radius
	mmcOpts += "-length=%d " % length
	mmcOpts += "-seed=%d " % seed
	
	if impl.lower() == 'mmc':
		from icecube import sim_services, c2j_icetray, mmc_icetray
		jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
		jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
		jvmOpts.append("-Xms256m")
		jvmOpts.append("-Xmx512m")
		jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")
		
		jvm = c2j_icetray.I3JavaVM(jvmOpts)
		return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)
	else:
		from icecube import sim_services, PROPOSAL_icetray
		return PROPOSAL_icetray.I3PropagatorServicePROPOSAL(mmcOpts)

@traysegment
def GenerateBundles(tray, name, Generator=None, Propagator=None,
    RunNumber=1, NEvents=100, RandomService=None,
    GCDFile='/data/sim/sim-new/downloads/GCD_31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz',
    FromTime=dataclasses.I3Time(55380),
    ToTime=dataclasses.I3Time(55380)):
	"""
	Generate muon bundles from a parametrization.
	
	:param Generator: an instance of I3MuonGun.Generator
	:param Propagator: an instance of I3
	"""
	
	from icecube import icetray, dataclasses
	from icecube import sim_services, MuonGun
	
	tray.AddModule("I3InfiniteSource",name+"_streams",
	    Prefix=GCDFile, Stream=icetray.I3Frame.DAQ)
	
	tray.AddModule("I3MCEventHeaderGenerator",name+"_gen_header",
	    Year=FromTime.utc_year, DAQTime=FromTime.utc_daq_time,
	    RunNumber=RunNumber, EventID=1, IncrementEventID=True)
	
	# modify the header if necessary by generating a random time
	# between StartTime and EndTime
	def GenRandomTime(frame, StartTime, EndTime, RandomService):
		header = frame["I3EventHeader"]
		del frame["I3EventHeader"]

		mjdStart = StartTime.mod_julian_day_double
		mjdEnd	 = EndTime.mod_julian_day_double
		
		mjd = mjdStart + (RandomService.uniform(0.,1.)*(mjdEnd-mjdStart))
		eventTime = dataclasses.I3Time()
		eventTime.set_mod_julian_time_double(mjd)
		header.start_time = eventTime
		frame["I3EventHeader"] = header
		
		if "DrivingTime" in frame:
			del frame["DrivingTime"]
			frame["DrivingTime"] = dataclasses.I3Time(eventTime)
		
	if FromTime != ToTime:
		if RandomService is None:
			raise ValueError("You must specify a random service!")
		tray.AddModule(GenRandomTime, name+"_GenRandomTime",
			StartTime=FromTime,
			EndTime=ToTime,
			RandomService=RandomService,
			Streams=[icetray.I3Frame.DAQ])
	
	tray.AddModule('I3MuonGun::GeneratorModule',name+"_gen_bundles",
	    Generator=NEvents*Generator, Propagator=Propagator)
	# tray.AddModule('Dump', name+'dump')
	
