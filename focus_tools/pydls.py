import os
import sys
import shutil

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

def dls76(args=[]):
    spgr_dat  = os.path.join(drc, "..", "resources", "spgr.dat")

    dls76_exe = '_dls76.x'

    if not os.path.exists(spgr_dat):
        print("Cannot find", spgr_dat)
        sys.exit()

    clean()

    try:
        inp = args[0]       
    except IndexError:
        inp = "dls76.inp"

    try:
        out = args[1]       
    except IndexError:
        out = "dls76.out"

    if not os.path.exists(inp):
        print("Cannot find", inp)
        sys.exit()

    inp = os.path.abspath(inp)
    out = os.path.abspath(out)

    shutil.copyfile(inp,      "fort.7")
    shutil.copyfile(spgr_dat, "fort.10")

    sp.call([dls76_exe,])

    move_files()
    clean()

def dls76_entry():
    args = sys.argv[1:]
    dls76(args)

if __name__ == '__main__':
    dls76_entry_point()