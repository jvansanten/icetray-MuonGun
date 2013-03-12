#!/usr/bin/env python

from icecube import icetray, dataio, MuonGun

generator = MuonGun.StaticSurfaceInjector()

frame = icetray.I3Frame()
frame['Generator'] = generator

f = dataio.I3File('muongun_serialization_test.i3', 'w')
f.push(frame)
f.close()
