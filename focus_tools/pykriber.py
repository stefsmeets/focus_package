import os
import sys
import shutil
from pathlib import Path
import subprocess as sp
from subprocess import PIPE, STDOUT
import platform

if sys.version_info < (3, 10):
    from importlib_resources import files
else:
    from importlib.resources import files


def clean():
    filenames = ["symdat", "distdat", "coseq"]
    for filename in filenames:
        if Path(filename).exists():
            os.remove(filename)


def setup():
    filenames = ["symdat", "distdat", "coseq"]

    resources = files('focus_tools.resources')
    
    for filename in filenames:
        src = resources / filename
        target = Path('.').absolute() / filename

        shutil.copyfile(src, target)


def prepare():
    clean()
    setup()

    filenames = ["symdat", "distdat", "coseq"]
    for filename in filenames:
        assert Path(filename).exists()

    if not Path("strudat").exists():
        print(">> Warning: Cannot find 'strudat' file\n")

    sysname = platform.system()

    if sysname == 'Linux':
        bin_dir = files('focus_tools.bin_linux')
    elif sysname == 'Darwin':
        bin_dir = files('focus_tools.bin_osx')
    elif sysname == 'Windows':
        bin_dir = files('focus_tools.bin_windows')
    else:
        raise RuntimeError(f'Unknown platform: {sysname}')

    return bin_dir / '_kriber.x'


def move(fname, target):
    if Path(target).exists():
        os.remove(target)
    os.rename(fname, target)


def extract_all_keys_from_strudat():
    keys = []
    try:
        strudat = open("strudat")
    except OSError:
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
        out = p.communicate(input=key.encode())
        stdout, stderr = out

        if stdout:
            stdout = stdout.decode()

            for line in stdout.split('\n'):
                if line.startswith('Traceback'):
                    print(stdout)
                    exit()
                elif "ERROR" in line:
                    raise RuntimeError(f"{key}  ->  KRIBER {line}")
        if stderr:
            print(stderr.decode())

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
        out = p.communicate(input=key.encode())
        stdout, stderr = out

        if stdout:
            stdout = stdout.decode()

            for line in stdout.split('\n'):
                if line.startswith('Traceback'):
                    print(stdout)
                    exit()
                elif "ERROR" in line:
                    raise RuntimeError(f"{key}  ->  KRIBER {line}")
        if stderr:
            print(stderr.decode())

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