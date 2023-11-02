![PyPI - Python Version](https://img.shields.io/pypi/pyversions/focus-package)
![PyPI](https://img.shields.io/pypi/v/focus-package.svg?style=flat)


FOCUS
=====

FOCUS is a program for structure determination of zeolites. It combines automatic Fourier recycling with a specialized topology search specific to zeolites or zeolite-like materials. This makes it ideal for structure determination from powder or electron diffraction data, because the framework search can make up for the loss of information from peak overlap or dynamical effects. This package also includes KRIBER, DLS76, and several python helper scripts that work together to analyze the results, optimize the framework geometry, and produce CIF files.

Installation
============

If you are a python user, you can install via pip:

    pip install focus-package

Alternatively, you can download one of the packages from the [releases](https://github.com/stefsmeets/focus_package/releases/latest) page.


Programs included
=================

 - focus
 - kriber
 - sginfo
 - dls76       
 - kriber      
 - fo2cif      
 - fo2hist     
 - fo2strudat  
 - strudat2cif 
 - dlsall      
 - cdlsall     
 - multifocal  

Please see the [manual](../master/manuals) for instructions.

Installation from source
========================

Extra dependencies for installation from source:

Build scripts:
 - [./focus/focus_build_script.sh](../master/focus/focus_build_script.sh)
 - [./src/kriber/makefile](../master/src/kriber/makefile)
 - [./src/dls76/makefile](../master/src/dls76/makefile)
 
References
==========

- S. Smeets, L. B. McCusker, C. Baerlocher, E. Mugnaioli, and U. Kolb. [Using FOCUS to solve zeolite structures from three-dimensional electron diffraction data](http://dx.doi.org/10.1107/S0021889813014817). J. Appl. Crystallogr., 46:1017–1023, 2013

- R. W. Grosse-Kunstleve, L. B. McCusker, C. Baerlocher. [Zeolite structure determination from powder diffraction data: applications of the FOCUS method](http://dx.doi.org/10.1107/S0021889899003453). J. Appl. Crystallogr., 32:536-542, 1999

- R. W. Grosse-Kunstleve, L. B. McCusker, C. Baerlocher. [Powder diffraction data and crystal chemical information combined in an automated structure determination procedure for zeolites](http://dx.doi.org/10.1107/S0021889897005013). J. Appl. Crystallogr., 30:985-995, 1997
