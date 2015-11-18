#!/usr/bin/env python

import os
import sys

import subprocess as sp

drc = os.path.abspath(os.path.dirname( __file__ )) # get path of program

filemap = {
	"fort.7"  : "dls76.inp",
	"fort.8"  : "nfilea.inp",
	"fort.10" : "spgr.dat",
	"fort.11" : "powd10.inp",
	"fort.12" : "dls76.out" }

def clean():
	fort = ["fort.7", "fort.8", "fort.10", "fort.11", "fort.12"]
	for f in fort:
		if os.path.exists(f):
			os.remove(f)

def move_files():
	fort = ["fort.8", "fort.11", "fort.12"]
	for f in fort:
		target = filemap[f]
		if os.path.exists(target):
			os.remove(target)
		if os.path.exists(f):
			os.rename(f, target)

def main():
	dls76_exe = os.path.join(drc, "dls_f", "dls76.x")
	spgr_dat  = os.path.join(drc, "dls_f", "spgr.dat")

	assert os.path.exists(dls76_exe)
	assert os.path.exists(spgr_dat)
	
	clean()

	args = sys.argv[1:]

	try:
		inp = args[0]		
	except IndexError:
		inp = "dls76.inp"

	try:
		out = args[1]		
	except IndexError:
		out = "dls76.out"

	assert os.path.exists(inp)

	inp = os.path.abspath(inp)
	out = os.path.abspath(out)

	# print inp
	# print out
	# print spgr_dat

	os.link(inp, "fort.7")
	os.link(out, "fort.10")

	sp.call([dls76_exe,])

	move_files()
	clean()

if __name__ == '__main__':
	main()