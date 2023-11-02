import os
import sys
import shutil
from pathlib import Path

import subprocess as sp

import platform
if sys.version_info < (3, 10):
    from importlib_resources import files
else:
    from importlib.resources import files


filemap = {
    "fort.7"  : "dls76.inp",
    "fort.8"  : "nfilea.inp",
    "fort.10" : "spgr.dat",
    "fort.11" : "powd10.inp",
    "fort.12" : "dls76.out" }

def clean():
    fort_filenames = ["fort.7", "fort.8", "fort.10", "fort.11", "fort.12"]
    for filename in fort_filenames:
        if Path(filename).exists():
            os.remove(filename)

def move_files():
    fort_filenames = ["fort.8", "fort.11", "fort.12"]
    for filename in fort_filenames:
        target = filemap[filename]
        if Path(target).exists():
            os.remove(target)
        if Path(filename).exists():
            os.rename(filename, target)

def dls76(args=[]):
    sysname = platform.system()

    if sysname == 'Linux':
        bin_dir = files('focus_tools.bin_linux')
    elif sysname == 'Darwin':
        bin_dir = files('focus_tools.bin_osx')
    elif sysname == 'Windows':
        bin_dir = files('focus_tools.bin_windows')
    else:
        raise RuntimeError(f'Unknown platform: {sysname}')

    dls76_exe = bin_dir / '_dls76.x'

    spgr_dat  = files('focus_tools.resources') / "spgr.dat"

    if not Path(spgr_dat).exists():
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

    if not Path(inp).exists():
        print("Cannot find", inp)
        sys.exit()

    inp = Path(inp).absolute()
    out = Path(out).absolute()

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