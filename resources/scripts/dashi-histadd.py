#!/usr/bin/env python

"""
Add all histograms in a group of files.
"""

import dashi
import tables
import os, sys

from optparse import OptionParser
parser = OptionParser(usage="%prog [OPTIONS] infiles outfile")
parser.add_option("--dest", dest="group", default="/", help="Copy datasets to this path in output file")
parser.add_option("--overwrite", dest="overwrite", default=False, action="store_true", help="Overwrite groups that exist in destination file")
parser.add_option("--norm", dest="norm", default=None, type=int, help="")
parser.add_option("--add", dest="add", default=False, action="store_true", help="Add histograms to those that already exist in destination file rather than overwritting them")

opts, args = parser.parse_args()
if len(args) < 2:
	parser.error("You must specify at least one output and one input file")
infiles, outfile = args[:-1], args[-1]

if opts.overwrite and opts.add:
	parser.error("--overwrite and --add are mutually exclusive!")

if os.path.exists(outfile) and not (opts.overwrite or opts.add):
	parser.error("%s already exists! Specify --overwrite or --add to modify it" % outfile)

# If adding, read in the current contents of the destination group
# before the remaining files
if opts.add:
	infiles = [outfile] + infiles
	root = opts.group
else:
	root = '/'

from collections import defaultdict
hists = dict()
counts = defaultdict(int)

for fname in infiles:
	print(fname)
	with tables.openFile(fname) as hdf:
		for group in hdf.walkNodes(where=root, classname='Group'):
			if 'ndim' in group._v_attrs: # a dashi histogram
				path = group._v_pathname
				h = dashi.histload(hdf, path)
				if root != '/' and path.startswith(root):
					path = path[len(root):]
				if path in hists:
					hists[path] += h
				else:
					hists[path] = h
				if 'count' in group._v_attrs:
					count = group._v_attrs['count']
				else:
					count = 1
				counts[path] += count
	# After the first file (possibly the destination), flip back to the root
	root = '/'

with tables.openFile(outfile, 'a') as hdf:
	if opts.group != "/":
		try:
			hdf.getNode(opts.group)
			hdf.removeNode(opts.group, recursive=True)
			print("removed %s" % (opts.group))
		except:
			pass
		hdf.createGroup("/", opts.group[1:])
	for path, h in hists.items():
		where = os.path.dirname(path)
		name = os.path.basename(path)
		if len(where) > 0:
			try:
				hdf.getNode(opts.group+where)
			except tables.NoSuchNodeError:
				hdf.createGroup(opts.group, where[1:], createparents=True)
		else:
			where = '/'
		dashi.histsave(h, hdf, opts.group+where, name)
		if opts.norm is not None:
			counts[path] = opts.norm
		hdf.getNode(opts.group+where+"/"+name)._v_attrs['count'] = counts[path]
		print("%s: %d entries" % (opts.group+where+"/"+name, counts[path]))

