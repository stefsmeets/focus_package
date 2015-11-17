#!/usr/bin/env sh

scons_ver="2.2.0"
scons="scons-$scons_ver"
fftw="fftw-3.3.2"

# download and extract focus
#curl -o focus-src.tar.gz http://cci.lbl.gov/~rwgk/focus/focus-src.tar.gz
#tar -xvf focus-src.tar.gz

# download libtbx
svn co https://cctbx.svn.sourceforge.net/svnroot/cctbx/trunk/libtbx libtbx --non-interactive --trust-server-cert

# download and install scons
curl -Lo $scons.tar.gz http://sourceforge.net/projects/scons/files/scons/$scons_ver/$scons.tar.gz/download
tar -xvf $scons.tar.gz
rm -rf $scons.tar.gz

mv $scons scons

mkdir build
cd build

# download and build fftw
prefix=`pwd`

curl -o $fftw.tar.gz http://www.fftw.org/$fftw.tar.gz
tar -xvf $fftw.tar.gz

cd $fftw
./configure --prefix="$prefix/base" --enable-single

make $*
make install

cd $prefix
rm -rf $fftw
rm -rf $fftw.tar.gz

# build focus
python ../libtbx/configure.py focus
make

echo ""
echo "All done! To recompile, run 'make redo' in ./build/"
echo "Executables for sginfo and focus can be found in ./build/exe/"
echo ""

