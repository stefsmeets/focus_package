#!/usr/bin/env

from setuptools import setup, find_packages
from os import path

setup(
    name="focus_package",
    version="0.1.0",
    description="FOCUS package including DLS-76 and KRIBER",

    author="Stef Smeets",
    author_email="stef.smeets@mat.ethz.ch",
    license="GPL",
    url="https://github.com/stefsmeets/focus_package",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["focus_package"],

    install_requires=[],

    package_data={
        "": ["LICENCE", "readme.md"],
        "focus_package": ["*.py", "dls_f/*", "kriber_f/*"]
    },

    scripts=["focus_package/bin/cdlsall",
             "focus_package/bin/coseq_cmp",
             "focus_package/bin/coseq_reduce",
             "focus_package/bin/dlsall",
             "focus_package/bin/fo2hist",
             "focus_package/bin/focus",
             "focus_package/bin/focus2cif",
             "focus_package/bin/focus2cm",
             "focus_package/bin/focus2strudat",
             "focus_package/bin/focus_check",
             "focus_package/bin/genseq",
             "focus_package/bin/get_coseq",
             "focus_package/bin/section",
             "focus_package/bin/sginfo",
             "focus_package/bin/strudat2cif"],

    entry_points={
        'console_scripts': [
            'multifocal = focus_package.multifocal:main',
            'dls76      = focus_package.pydls:main',
            'kriber     = focus_package.pykriber:main',
        ]
    }
)
