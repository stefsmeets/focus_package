FOCUS
=====

FOCUS is a program for structure determination of zeolites. It combines automatic Fourier recycling with a specialized topology search specific to zeolites or zeolite-like materials. This makes it ideal for structure determination from powder of electron diffraction data, because the framework search can make up for the loss of information from peak overlap or dynamical effects. This package also includes KRIBER, DLS76, and several python helper scripts that work together to analyze the results, optimize the framework geometry, and produce CIF files.

Installation
============

    pip install focus-package

### Windows

The windows version includes all dependencies and can be used as is. In case Windows complains, it is sometimes needed to install [Visual C++ Redistributable for Visual Studio 2015](https://www.microsoft.com/en-us/download/details.aspx?id=48145) (vc_redist.x86.exe)

### OSX

The OSX version of DLS76 has been compiled using brew gfortran [homebrew](http://brew.sh/)  
To install gfortran:
    
    brew install gfortran

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

(please see the [manual](../master/manuals))

Installation from source
========================

Extra dependencies for installation from source:

### OS X

- command line tools (xcode-select install)
- gfortran

### Linux

- svn
- curl

Extract using:

    tar -xvzf focus_package-src.tar.gz

Build scripts:
 - [./focus/focus_build_script.sh](../master/focus/focus_build_script.sh)
 - [./src/kriber/makefile](../master/src/kriber/makefile)
 - [./src/dls76/makefile](../master/src/dls76/makefile)
 
References
==========
 
S. Smeets, L. B. McCusker, C. Baerlocher, E. Mugnaioli, and U. Kolb. [Using FOCUS to solve zeolite structures from three-dimensional electron diffraction data](http://dx.doi.org/10.1107/S0021889813014817). J. Appl. Crystallogr., 46:1017–1023, 2013

R. W. Grosse-Kunstleve, L. B. McCusker, C. Baerlocher. [Zeolite structure determination from powder diffraction data: applications of the FOCUS method](http://dx.doi.org/10.1107/S0021889899003453). J. Appl. Crystallogr., 32:536-542, 1999

R. W. Grosse-Kunstleve, L. B. McCusker, C. Baerlocher. [Powder diffraction data and crystal chemical information combined in an automated structure determination procedure for zeolites](http://dx.doi.org/10.1107/S0021889897005013). J. Appl. Crystallogr., 30:985-995, 1997
