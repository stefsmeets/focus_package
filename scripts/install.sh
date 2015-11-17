#!/usr/bin/env sh


base=`pwd`

kriber_dir="kriber_f"
dls_dir="dls_f"


mkdir $base/bin

cp build/exe/*         $base/bin/

cp focus/focus_check   $base/bin/
cp focus/multifocal.py $base/bin/multifocal

cp utils/coseq_cmp     $base/bin/
cp utils/coseq_reduce  $base/bin/
cp utils/section       $base/bin/
cp utils/genseq        $base/bin/

cp scripts/*           $base/bin/

ln kriber_f/kriber.csh $base/bin/kriber

ln dls_f/dls76.csh     $base/bin/dls76


echo ""
echo "================================================================================"
echo ""
echo " All done!!"
echo ""
echo " ** for bash/zsh, add these lines to your ~/.zshrc or ~/.bashrc file:"
echo ""
echo "        export PATH=\"\${PATH}:$base/bin\""
echo "        export COSEQ_DB=\"$base/$kriber_dir/coseq\""
echo
echo " ** for csh/tcsh, add these lines to your ~/.cshrc or ~/.tcshrc file:"
echo ""
echo "        set path=(\$path $base/bin)"
echo "        setenv COSEQ_DB $base/$kriber_dir/coseq"
echo ""