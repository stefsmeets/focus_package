#!/usr/bin/env bash

base=`pwd`

rm -rf ./bin/

cd $base/focus
sh focus_build_script.sh

cd $base/src/kriber
make clean
make

cd $base/src/dls76
make clean
make

# cd $base/utils
# make

cd $base

mkdir bin/
mkdir bin/bin

mv src/kriber/kriber.x    bin/bin
mv src/dls76/dls76.x      bin/bin
mv focus/build/exe/focus  bin/
mv focus/build/exe/sginfo bin/

echo
echo "To install: python setup.py install"
echo
