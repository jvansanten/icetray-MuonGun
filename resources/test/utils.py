
import random
from os.path import expandvars
def get_propagator(radius=800, length=1600, seed=random.randint(0, (1<<32) - 1)):
	mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
	mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
	mmcOpts += expandvars("-mediadef=$I3_BUILD/mmc-icetray/resources/mediadef ")
	mmcOpts += "-radius=%d " % radius
	mmcOpts += "-length=%d " % length
	mmcOpts += "-seed=%d " % seed
	try:
		from icecube import PROPOSAL_icetray
		return PROPOSAL_icetray.I3PropagatorServicePROPOSAL(mmcOpts)
	except ImportError:
		pass
	
	try:
		from icecube import c2j_icetray, mmc_icetray
		jvmOpts = icetray.vector_string()
		jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
		jvmOpts.append("-Xms256m")
		jvmOpts.append("-Xmx512m")
		jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")
		
		jvm = c2j_icetray.I3JavaVM(jvmOpts)
		return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)
	except ImportError:
		pass
	
	print('Neither PROPOSAL nor MMC found!')
	return None
