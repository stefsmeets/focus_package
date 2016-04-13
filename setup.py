#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path
import sys

if sys.platform == "win32":
    package_dir={"bin": "bin_windows/bin"}
elif sys.platform == "darwin":
    package_dir={"bin": "bin_osx/bin"}
elif sys.platform == "linux2":
    package_dir={"bin": "bin_linux/bin"}
else:
    raise RuntimeError

print package_dir

setup(
    name="focus_package",
    version="1.1",
    description="FOCUS package including DLS-76 and KRIBER",

    author="Stef Smeets",
    author_email="stef.smeets@mmk.su.se",
    maintainer="Stef Smeets",
    maintainer_email="stef.smeets@mmk.su.se",
    license="GPL",
    url="https://github.com/stefsmeets/focus_package",

    classifiers=[
        'Programming Language :: Python :: 2.7',
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
    }
)

if sys.platform == "win32":
    print """
 >> Installation complete. You can now remove this directory.

Put these files in a location on your system path:

    ./bin_windows/libfftw3f-3.dll
    ./bin_windows/focus.exe
    ./bin_windows/sginfo.exe
"""
else:    
    print """
 >> Installation complete. You can remove this directory.

Put these files in a location on your system path:
"""
    if sys.platform == "darwin"
        print """     
        ./bin_osx/focus
        ./bin_osx/sginfo"""
    elif sys.platform == "linux2":
        print """     
        ./bin_linux/focus
        ./bin_linux/sginfo"""

    



