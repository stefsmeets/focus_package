#!/usr/bin/env sh

scons_ver="2.4.1"
scons="scons-$scons_ver"
fftw="fftw-3.3.2"
# New versions of libtbx break compilation of FOCUS, because boost cannot be found
libtbxrev=20000

if [ ! -d libtbx ];
then
	# download libtbx
	# svn co https://cctbx.svn.sourceforge.net/svnroot/cctbx/trunk/libtbx libtbx --non-interactive --trust-server-cert -r$libtbxrev
    svn co https://github.com/cctbx/cctbx_project/trunk/libtbx libtbx --non-interactive --trust-server-cert -r$libtbxrev
fi

if [ ! -d scons ];
then
	# download and install scons
	curl -Lo $scons.tar.gz http://sourceforge.net/projects/scons/files/scons/$scons_ver/$scons.tar.gz/download
	tar -xvf $scons.tar.gz
	rm -rf $scons.tar.gz

	mv $scons scons
fi

mkdir build
cd build

if [ ! -f ./base/bin/fftwf-wisdom ];
then
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
fi

# build focus
python ../libtbx/configure.py src
make

echo ""
echo "All done! To recompile, run 'make redo' in ./build/"
echo "Executables for sginfo and focus can be found in ./build/exe/"
echo ""

