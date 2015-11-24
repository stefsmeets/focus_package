#!/usr/bin/env python

import os
import sys
import shutil

import subprocess as sp

drc = os.path.abspath(os.path.dirname( __file__ )) # get path of program

print "RAWR"

def clean():
	files = ["symdat", "distdat", "coseq"]
	for f in files:
		if os.path.exists(f):
			os.remove(f)

def setup():
	files = ["symdat", "distdat", "coseq"]
	for f in files:

		src = os.path.join(drc, "kriber_f", f)
		target = os.path.join(os.path.abspath("."), f)

		shutil.copyfile(src, target)

def main():
	kriber_exe = os.path.join(drc, "kriber_f", "kriber.x")

	clean()
	setup()

	assert os.path.exists(kriber_exe)
	files = ["symdat", "distdat", "coseq"]
	for f in files:
		assert os.path.exists(f)

	sp.call([kriber_exe] + sys.argv[1:])

	clean()

if __name__ == '__main__':
	main()