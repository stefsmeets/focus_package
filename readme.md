Quick installation
==================

Download the FOCUS package from https://www.crystal.mat.ethz.ch/Software/


Extract using

    tar -xvzf focus_package-linux.tar.gz

or

    tar -xvzf focus_package-osx.tar.gz

Install with python:

    python setup.py install


The OSX version of DLS76 have been compiled using brew gfortran (<http://brew.sh/>)
To install gfortran:
    
    brew install gfortran


After the install step, you can remove the whole directory. 
The scripts have been tested to work on recent versions of OS X and Ubuntu.


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

Build using:

    sh build_all.sh

Install with python:

    python setup.py install



