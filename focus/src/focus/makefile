objects = atominfo.o cputime.o evalfw.o findpnnb.o frame.o io.o\
          lattice.o lib.o matrix.o posnmx.o psearch.o ranmar.o recycle.o\
	  scatfact.o setaos.o sgextend.o signal.o sitefcal.o tetrahed.o\
          trialoop.o trpcmplx.o xtal.o nodsurf.o main.o

sources = atominfo.c cputime.c evalfw.c findpnnb.c frame.c io.c\
          lattice.c lib.c matrix.c posnmx.c psearch.c ranmar.c recycle.c\
	  scatfact.c setaos.c sgextend.c signal.c sitefcal.c tetrahed.c\
          trialoop.c trpcmplx.c xtal.c main.c fourier.c nodsurf.c

program = focus

CC = cc
CCFLAGS = -O
CCDEFS  =#-DNO_EXTRA_CHECKS#-DPROVOKE_SEG_FAULT#-DCOSEQ_TRACE#-DFP_USE_DOUBLE

# Paragon:
#CCFLAGS = -O -nx

LOCALINC = -I$(LOCAL)/include
LOCALLIB = -L$(LOCAL)/lib

CCMOREF = $(CCFLAGS) $(LOCALINC)

SHELL=/bin/sh

what:
	@echo 'make what?'

focus.nofft: $(objects) fourier.nofft.o
	$(CC) $(CCMOREF) $(objects) fourier.nofft.o -o $(program).nofft $(LOCALLIB) -lsginfo -lm

$(program).fftw: fourier.fftw.o $(objects)
	$(CC) $(CCMOREF) $(objects) fourier.fftw.o -o $(program).fftw $(LOCALLIB) -lsginfo -lrfftw -lfftw -lm

$(program).fftpack: fourier.fftpack.o $(objects)
	$(CC) $(CCMOREF) $(objects) fourier.fftpack.o -o $(program).fftpack $(LOCALLIB) -lsginfo -lrfftw -lfftw -lfftpack41 $(LIBFTN) -lm

$(program).fftnsgi: fourier.fftnsgi.o $(objects)
	$(CC) $(CCMOREF) $(objects) fourier.fftnsgi.o -o $(program).fftnsgi $(LOCALLIB) -lsginfo -lrfftw -lfftw -lfftpack41 -lcomplib.sgimath -lftn -lm

$(program).decunix: fourier.decdxml.o $(objects)
	$(CC) $(CCMOREF) $(objects) fourier.decdxml.o -o $(program).decunix $(LOCALLIB) -lsginfo -lrfftw -lfftw -lfftpack41 -ldxml -lm

$(program).hpux: fourier.fftpack.o $(objects)
	$(CC) $(CCMOREF) $(objects) fourier.fftpack.o -o $(program).hpux    $(LOCALLIB) -lsginfo -lrfftw -lfftw -lfftpack41 -lm -lcl

fourier.nofft.o: fourier.c
	$(CC) $(CCMOREF) $(CCDEFS) -o fourier.nofft.o -c fourier.c

fourier.fftw.o: fourier.c
	$(CC) $(CCMOREF) $(CCDEFS) -o fourier.fftw.o -c fourier.c -DFFTfftw

fourier.fftpack.o: fourier.c
	$(CC) $(CCMOREF) $(CCDEFS) -o fourier.fftpack.o -c fourier.c -DFFTfftpack -DFFTfftw

fourier.fftnsgi.o: fourier.c
	$(CC) $(CCMOREF) $(CCDEFS) -o fourier.fftnsgi.o -c fourier.c -DFFTfftpack -DFFTfftw -DFFTsgimath

fourier.decdxml.o: fourier.c
	$(CC) $(CCMOREF) $(CCDEFS) -o fourier.decdxml.o -c fourier.c -DFFTfftpack -DFFTfftw -DFFTdecdxml

.c.o:
	$(CC) $(CCMOREF) $(CCDEFS) -c $*.c
        
lint:
	lint -DFFTfftpack -DFFTsgimath $(sources) $(LOCALINC)

focusall.c: $(sources)
	echo '#define A_AppGlobal__' > focusall.c
	cat $(sources) >> focusall.c

atominfo: atominfo.c
	$(CC) $(CCMOREF) -o atominfo atominfo.c -DSTAND_ALONE -lm

lattice: lattice.c
	$(CC) $(CCMOREF) -o lattice lattice.c -DSTAND_ALONE -lm

cputime: cputime.c
	$(CC) $(CCMOREF) -o cputime cputime.c -DSTAND_ALONE -lm

matrix: matrix.c
	$(CC) $(CCMOREF) -o matrix matrix.c -DSTAND_ALONE -lm

ranmar: ranmar.c
	$(CC) $(CCMOREF) -o ranmar ranmar.c -DSTAND_ALONE_1 -lm


clean:
	rm -f $(objects) focusall.c
	rm -f fourier.nofft.o fourier.fftw.o fourier.fftpack.o
	rm -f fourier.fftnsgi.o fourier.decdxml.o
	rm -f focus.fftpack focus.fftnsgi focus.decunix focus.hpux core
	rm -f atominfo lattice cputime matrix ranmar
