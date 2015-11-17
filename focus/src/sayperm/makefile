#CC = cc
#CCFLAGS = -O2 -woff 1552
#HOMELIB = -L$(HOME)/focus-src/local/sgi/lib
CC = gcc
CCFLAGS = -O2
HOMELIB = -L../../libtbx
objects =  sayperm.o golay.o sayre.o subrout.o atomfac.o
sayperm:$(objects)
	$(CC) $(CCFLAGS) -o sayperm $(objects) $(HOMELIB) -lsginfo -lm 
.c.o:
	$(CC) $(CCFLAGS) -c $*.c
golay:
	$(CC) $(CCFLAGS) golay.c -o golay 
clean:
	rm -f $(objects) sayperm.o sayperm golay
