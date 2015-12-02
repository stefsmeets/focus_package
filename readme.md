Quick installation
==================

Compiled versions and source distribution are available from https://github.com/stefsmeets/focus_package/releases


Extract using

    tar -xvzf focus_package-linux.tar.gz

or

    tar -xvzf focus_package-osx.tar.gz

Install with python:

    python setup.py install

This will install all files into a directory on your system path. After this, you can remove the installation directory. 

The OSX version of DLS76 has been compiled using brew gfortran (<http://brew.sh/>)
To install gfortran:
    
    brew install gfortran

These scripts have been tested to work on recent versions of OS X and Ubuntu.


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



