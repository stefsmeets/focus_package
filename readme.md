Quick installation
==================

Download the FOCUS package from the [releases page](https://github.com/stefsmeets/focus_package/releases).

1.	Unzip or extract using:  `tar -xvzf focus-osx.tar.gz`
2.	Install with python: `python setup.py install`

Focus and all included scripts have been tested to work on recent versions of Windows, OS X, Ubuntu.

### Windows

Windows users need to install [Python 2.7](https://www.python.org/downloads/)  
The Windows version also needs [Visual C++ Redistributable for Visual Studio 2015](https://www.microsoft.com/en-us/download/details.aspx?id=48145) (vc_redist.x86.exe)

### OSX

The OSX version of DLS76 has been compiled using brew gfortran [homebrew](http://brew.sh/)  
To install gfortran:
    
    brew install gfortran

Tools included
==============

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

(please see the [manual](../blob/master/manuals))

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
 - [./focus/focus_build_script.sh](../blob/master/focus/focus_build_script.sh)
 - [./src/kriber/makefile](../blob/master/src/kriber/makefile)
 - [./src/dls76/makefile](../blob/master/src/dls76/makefile)
