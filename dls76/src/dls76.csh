#! /bin/csh -f

set progn = $0
set progn = $progn:t

set DLSDIR = ~/bin/DLS76-02

set SPGR_DAT = $DLSDIR/spgr.dat
set DLSX = $DLSDIR/dls76.x

# If you never want powd10.inp you may change the Fortran code (see there)

# Chanel  I/O   File
#    7     I    dls76.inp
#    8     O    nfilea.inp
#   10     I    spgr.dat
#   11     O    powd10.inp
#   12     O    dls76.out

rm -f fort.7 fort.8 fort.10 fort.11 fort.12

if ($#argv > 0) then
  set dls76_inp = $1
else
  set dls76_inp = dls76.inp
endif

if ($#argv > 1) then
  set dls76_out = $2
else
  set dls76_out = dls76.out
endif

if (! -e $dls76_inp) then
  echo $progn': Error: Missing '$dls76_inp'\07'
  exit 1
endif
if (! -e $SPGR_DAT) then
  echo $progn': Error: Missing '$SPGR_DAT'\07'
  exit 1
endif

onintr CLEAN

ln -s $dls76_inp fort.7
ln -s $SPGR_DAT fort.10

$DLSX

CLEAN:

if (-e fort.8) mv fort.8 nfilea.inp
if (-e fort.11) mv fort.11 powd10.inp
if (-e fort.12) mv fort.12 $dls76_out

rm -f fort.7 fort.8 fort.10 fort.11 fort.12

# SGI only
if (`ls | grep 'tmp\.Faaa.*' | wc -l` != 0) then
  rm -f tmp.Faaa*
endif

exit 0
