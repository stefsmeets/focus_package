#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import shutil

import subprocess as sp
from subprocess import PIPE, STDOUT

drc = os.path.abspath(os.path.dirname( __file__ )) # get path of program


def clean():
    files = ["symdat", "distdat", "coseq"]
    for f in files:
        if os.path.exists(f):
            os.remove(f)


def setup():
    files = ["symdat", "distdat", "coseq"]
    for f in files:

        src = os.path.join(drc, "..", "resources", f)
        target = os.path.join(os.path.abspath("."), f)

        shutil.copyfile(src, target)


def prepare():
    kriber_exe = os.path.join(drc, "..", "bin", "kriber.x")
    kriber_exe_dev = os.path.join(drc, "..", "bin_windows", "bin", "kriber.x") # check developer path

    clean()
    setup()

    if not os.path.exists(kriber_exe):
        if not os.path.exists(kriber_exe_dev):
            print("Cannot find", kriber_exe)
            sys.exit()
        kriber_exe = kriber_exe_dev

    files = ["symdat", "distdat", "coseq"]
    for f in files:
        assert os.path.exists(f)

    if not os.path.exists("strudat"):
        print(">> Warning: Cannot find 'strudat' file\n")

    return kriber_exe


def move(fname, target):
    if os.path.exists(target):
        os.remove(target)
    os.rename(fname, target)


def extract_all_keys_from_strudat():
    keys = []
    try:
        strudat = open("strudat", "r")
    except IOError:
        print("Cannot find 'strudat' file")
        sys.exit()

    for line in strudat:
        if line.startswith("*"):
            keys.append(line.strip()[1:])
    return keys


def strudat2cif(args=[], keys=[], rename=True, verbose=True):
    from subprocess import PIPE, STDOUT
    
    if isinstance(keys, str):
        keys = keys.split()

    kriber_exe = prepare()
    
    if not keys:
        keys = extract_all_keys_from_strudat()
    
    for key in keys:
        p = sp.Popen(['kriber', 'reacs'] + args + ['wricif', 'exit'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        out = p.communicate(input=key)

        if out[0]:
            for line in out[0].split("\n"):
                if "ERROR" in line:
                    raise RuntimeError("{}  ->  KRIBER {}".format(key, line))
        if out[1]:
            print(out[1])

        move("structure.cif", key+".cif")
        if verbose:
            print(" >> Wrote file {}".format(key+".cif"))

    clean()


def strudat2dls(args=[], keys=[], verbose=True):
    from subprocess import PIPE, STDOUT
    
    if isinstance(keys, str):
        keys = keys.split()

    kriber_exe = prepare()
    
    if not keys:
        keys = extract_all_keys_from_strudat()

    for key in keys:
        p = sp.Popen(['kriber', 'reacs'] + args + ['wriid', 'exit'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        out = p.communicate(input=key)

        if out[0]:
            for line in out[0].split("\n"):
                if "ERROR" in line:
                    raise RuntimeError("{}  ->  KRIBER {}".format(key, line))
        if out[1]:
            print(out[1])

        # move("structure.cif", key+".cif")
        # print " >> Wrote file {}".format(key+".cif")

    clean()


def kriber(args=[]):
    kriber_exe = prepare()
    sp.call([kriber_exe] + args)
    clean()


def strudat2cif_entry():
    fns = sys.argv[1:]
    strudat2cif(args=fns)


def kriber_entry():
    fns = sys.argv[1:]
    kriber(args=fns) 


if __name__ == '__main__':
    run_kriber()
    # strudat2cif()