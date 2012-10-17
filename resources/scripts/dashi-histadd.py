#!/usr/bin/env python

"""
Add all histograms in a group of files.
"""

import dashi
import tables
import os, sys

infiles, outfile = sys.argv[1:-1], sys.argv[-1]

if os.path.exists(outfile):
	raise ValueError("%s already exists!" % outfile)

hists = dict()

for fname in infiles:
	with tables.openFile(fname) as hdf:
		for group in hdf.walkNodes(classname='Group'):
			if 'ndim' in group._v_attrs: # a dashi histogram
				path = group._v_pathname
				h = dashi.histload(hdf, path)
				if path in hists:
					hists[path] += h
				else:
					hists[path] = h

with tables.openFile(outfile, 'w') as hdf:
	for path, h in hists.iteritems():
		where = os.path.dirname(path)
		name = os.path.basename(path)
		if len(where) > 0:
			try:
				hdf.getNode(where)
			except tables.NoSuchNodeError:
				hdf.createGroup('/', where[1:], createparents=True)
		else:
			where = '/'
		dashi.histsave(h, hdf, where, name)

