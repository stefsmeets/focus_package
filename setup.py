#!/usr/bin/env python

from setuptools import setup, find_packages
import os, sys

if "win" in sys.argv:
    platform = "win32"
elif "macosx" in sys.argv:
    platform = "darwin"
elif "linux" in sys.argv:
    platform = "linux"
else:
    platform = sys.platform

if platform == "win32":
    package_dir={"bin": "bin_windows/bin"}
    scripts = [os.path.join("bin_windows", fn) for fn in  ("libfftw3f-3.dll", "focus.exe", "sginfo.exe")]
elif platform == "darwin":
    package_dir={"bin": "bin_osx/bin"}
    scripts = [os.path.join("bin_osx", fn) for fn in  ("focus", "sginfo")]
elif platform == "linux":
    package_dir={"bin": "bin_linux/bin"}
    scripts = [os.path.join("bin_linux", fn) for fn in  ("focus", "sginfo")]
else:
    raise RuntimeError

setup(
    name="focus_package",
    version="3.0.0",
    description="FOCUS package including DLS-76 and KRIBER",
    description_file = "README.md",

    author="Stef Smeets",
    author_email="s.smeets@esciencecenter.nl",
    license="GPLv2",
    url="https://github.com/stefsmeets/focus_package",

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering',
    ],

    package_dir=package_dir,

    packages=["focus_tools", "bin", "resources"],

    install_requires=[],

    package_data={
        "": ["LICENCE", "readme.md"],
        "focus_tools": ["*.py"],
        "resources": ["*"],
        "bin": ["*"]
    },

    entry_points={
        'console_scripts': [
            'multifocal  = focus_tools.multifocal:main',
            'dls76       = focus_tools.pydls:dls76_entry',
            'kriber      = focus_tools.pykriber:kriber_entry',
            'strudat2cif = focus_tools.pykriber:strudat2cif_entry',
            'fo2hist     = focus_tools.focus_tools:fo2hist_entry',
            'fo2cif      = focus_tools.focus_tools:fo2cif_entry',
            'fo2strudat  = focus_tools.focus_tools:fo2strudat_entry',
            'dlsall      = focus_tools.focus_tools:dlsall_entry',
            'cdlsall     = focus_tools.focus_tools:cdlsall_entry',
        ]
    },

    data_files = [
        ('bin', scripts)
    ],
)
