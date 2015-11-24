#!/usr/bin/env sh

base=`pwd`
kriber_dir="kriber_f"

# Installation and setup of kriber

cd src
make

cd $base
mkdir $kriber_dir

cp src/kriber.x $kriber_dir
cp resources/*  $kriber_dir

cd $kriber_dir

echo '#!/bin/csh

#Header of kriber.pas / kriber.c
#
#******************************************************************************
#*                                                                            *
#*   Program  KRIBER              Version 1.2       (February 2013)           *
#*   ===============                                                          *
#*                                                                            *
#*   Author:   Roland Bialek, Institut fuer Kristallografie und Petrografie   *
#*             ETH-Zentrum, CH-8092 Zuerich, Switzerland                      *
#*                                                                            *
#*   an interactive PASCAL program to                                         *
#*    - calculate distances and angles                                        *
#*    - generate input files for the programs DLS-76 and LOADAT (XRS-82)      *
#*    - calculate coordination sequences and loop configurations              *
#*                                                                            *
#*   Program Installation and Modification                                    *
#*   -------------------------------------                                    *
#*  The program is written in standard PASCAL and has been implemented and    *
#*  fully tested on a CYBER 855 running under NOS/VE 1.5.2. It has also been  *
#*  used successfully on a MicroVax II machine under VMS 5.2                  *
#*                                                                            *
#*  There is only one extension to the standard PASCAL used: The identifiers  *
#*  contain the underscore (_) and the dollar sign ($). If your PASCAL        *
#*  does not have this extension, you can remove these two signs with your    *
#*  editor.                                                                   *
#*                                                                            *
#*  Input:  The program runs interactively.                                   *
#*                                                                            *
#*  Output: The output is written to the terminal screen (longest line is     *
#*          80 columns wide).                                                 *
#*                                                                            *
#*  Files:  SYMDAT, STRUDAT, DISTDAT and COSEQ are input files, DLSINPUT,     *
#*          LOADATINPUT, LIST and STRUDATNEW are output files.                *
#*                                                                            *
#*  Arrays: Following constants can be changed in the constant declaration    *
#*          part:                                                             *
#*          maxanzahlatomlagen    :   maximum number of atom positions        *
#*          maxanzahlbindungen    :   maximum number of bonds per atom        *
#*          maxanzahlsymeqkarten  :   maximum number of symeq cards           *
#*          maxindepwinkel        :   maximum number of independent angles    *
#*          maxindepdistanzen     :   maximum number of independent distances *
#*          maxanzahlgegwinkel    :   maximum number of prescribed angles     *
#*          maxanzahlgegdistanzen :   maximum number of prescribed distances  *
#*          maxanzahldistanzen    :   maximum number of calculated distances  *
#*                                                                            *
#******************************************************************************

umask 022

onintr Interrupt

# Where symdat, strudat, distdat, coseq and kriber.x are' > kriber.csh

echo "set xdir = $base/$kriber_dir" >> kriber.csh
echo "set KRIBERX = $base/$kriber_dir/kriber.x" >> kriber.csh

echo '
# Executable
# set KRIBERX = ./kriber.x

# kriber Input files
set kif = (symdat strudat distdat coseq)
# "ln -s" list
set lnl = (0      0       0       0    )

set i=0
foreach f ($kif)
  @ i++
  if (! -e $f) then
    ln -s $xdir/$f $f
    set lnl[$i]=1
    echo "kriber: "$f" = "$xdir/$f
  else
    set lnl[$i]=0
    echo "kriber: "$f" = ./"$f
  endif
end

echo "Running "$KRIBERX" .."
$KRIBERX $*
set Xstatus=$status

UNLNS:

set i=0
foreach f ($kif)
  @ i++
  if ($lnl[$i] == 1) then
    rm -f $f
    set lnl[$i]=0
  endif
end

if (-e list) then
  set wclist=`wc -c list`
  if ("$wclist" =~ "0 list") rm -f list
endif

exit($Xstatus)


Interrupt:
  echo "Interrupt signal received."
  set Xstatus=1
  goto UNLNS
' >> kriber.csh

chmod +x kriber.csh
cd $base

unlink kriber
ln -s $base/$kriber_dir/kriber.csh kriber
