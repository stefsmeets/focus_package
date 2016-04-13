###########################################################################
# OSX   using gfortran                                                                 #
###########################################################################
CMPLR  = gfortran
LINKR  = gfortran
CFLAGS = -O2
LFLAGS = -O2 -static
# LFLAGS = -O2
LLIBS =   
###########################################################################
# OSX                                                                     #
###########################################################################
#CMPLR  = g77
#LINKR  = g77
#CFLAGS = -O2
#LFLAGS = -O2 -static
#LFLAGS = -O2
#LLIBS =   
###########################################################################
# OSX  intel Fortran                                                                   #
###########################################################################
#CMPLR  = ifort
#LINKR  = ifort
#CFLAGS = -O0 -mp
##LFLAGS = -O2 -static
#LFLAGS = -O0
#LLIBS =   

###########################################################################
# LINUX                                                                   #
###########################################################################
#CMPLR  = cc
#LINKR  = g77
#CFLAGS = -O2 -I../../lib_linux/include/
#LFLAGS = -O2 -static -L../../lib_linux/lib/ 
#LLIBS = -lreadline -lhistory -lncurses   


###########################################################################
# Don«t change anything below                                             #
###########################################################################

OBJS   =

.c.o :
	$(CMPLR) $(CFLAGS) -c $<

.f.o :
	$(CMPLR) $(CFLAGS) -c $<

dls76.x :  dls76.o
	$(LINKR) $(LFLAGS) -o $@ dls76.o $(OBJS) -lm $(LLIBS)

clean :
	rm -f *.o
	rm -f core
	rm -f dls76.x
	
