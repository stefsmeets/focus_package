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

mv kriber/src/kriber.x    bin/_kriber.x
mv dls76/src/dls76.x      bin/_dls76.x
mv focus/build/exe/focus  bin/
mv focus/build/exe/sginfo bin/

echo
echo "To install: python setup.py install"
echo
