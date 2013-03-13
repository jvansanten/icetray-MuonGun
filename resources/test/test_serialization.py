#!/usr/bin/env python

from icecube import icetray, dataio, MuonGun, phys_services

fname = 'muongun_serialization_test.i3'
generator = MuonGun.StaticSurfaceInjector()

frame = icetray.I3Frame()
frame['Generator'] = generator

f = dataio.I3File(fname, 'w')
f.push(frame)
f.close()

f = dataio.I3File(fname)
frame = f.pop_frame()
f.close()

newgenerator = frame['Generator']

