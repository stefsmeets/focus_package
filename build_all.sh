#!/usr/bin/env bash

base=`pwd`

cd $base/focus
sh focus_build_script.sh

cd $base/kriber_f
sh build_kriber_script.sh

cd $base/dls_f
sh build_dls_script.sh

cd $base/utils
make

cd $base

pkg="focus_package"
mkdir $pkg
mkdir $pkg/bin

mv kriber/kriber_f $pkg
mv dls76/dls_f $pkg

cp kriber/pykriber.py      $pkg/
cp dls76/pydls.py          $pkg/
cp utils/coseq_cmp         $pkg/bin/
cp utils/coseq_reduce      $pkg/bin/
cp utils/genseq            $pkg/bin/
cp utils/section           $pkg/bin/
cp scripts/*               $pkg/bin/
mv $pkg/bin/multifocal.py  $pkg/
mv focus/build/exe/*       $pkg/bin/
touch $pkg/__init__.py

echo
echo "To install: python setup.py install"
echo
