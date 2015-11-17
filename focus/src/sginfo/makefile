CC = cc
CCFLAGS = -O

SGINFOFLAGS =#-DNo_ListRotMxInfo

sginfo: sginfo.o sgvers.o sgclib.o sgio.o sgfind.o sghkl.o sgss.o sgharker.o sgutil.o sgssvfy.o xtal32ss.o sgcond.o sgasciio.o sghklvfy.o sgcharmx.o sggrid.o sgasym.o cnslib.o sgnorm.o
	$(CC) $(CCFLAGS) -o sginfo sginfo.o sgvers.o sgclib.o sgio.o sgfind.o sghkl.o sgss.o sgharker.o sgutil.o sgssvfy.o xtal32ss.o sgcond.o sgasciio.o sghklvfy.o sgcharmx.o sggrid.o sgasym.o cnslib.o sgnorm.o -lm

sginfo.o: sginfo.c sginfo.h
	$(CC) $(CCFLAGS) -c sginfo.c $(SGINFOFLAGS)

xtal32ss.o: xtal32ss.c sginfo.h
	$(CC) $(CCFLAGS) -c $*.c -DDEBUG_xtal32ss

sgtest: sgtest.o sgclib.o sgio.o sgfind.o sgutil.o sghkl.o sgcond.o sghklvfy.o sgcharmx.o
	$(CC) $(CCFLAGS) -o sgtest sgtest.o sgclib.o sgio.o sgfind.o sgutil.o sghkl.o sgcond.o sghklvfy.o sgcharmx.o

sgtest.o: sginfo.h

sgquick: sgquick.o sgclib.o sgio.o sgutil.o
	$(CC) $(CCFLAGS) -o sgquick sgquick.o sgclib.o sgio.o sgutil.o sgcharmx.o

sgquick: sginfo.h

sgdev: sgdev.o sgclib.o sgutil.o sgcharmx.o sgio.o
	$(CC) $(CCFLAGS) -o sgdev sgdev.o sgclib.o sgutil.o sgcharmx.o sgio.o

sgdev.o: sginfo.h


sgclib.o:   sginfo.h
sgio.o:     sginfo.h
sgfind.o:   sginfo.h
sghkl.o:    sginfo.h
sgss.o:     sginfo.h
sgharker.o: sginfo.h
sgutil.o:   sginfo.h
sgssvfy.o:  sginfo.h
xtal32ss.o: sginfo.h
sgcond.o:   sginfo.h
sgasciio.o: sginfo.h
sghklvfy.o: sginfo.h
sgcharmx.o: sginfo.h
sggrid.o:   sginfo.h
sgasym.o:   sginfo.h
cnslib.o:   sginfo.h
sgnorm.o:   sginfo.h

libsginfo: sgvers.o sgclib.o sgio.o sgfind.o sghkl.o sgss.o sgharker.o sgutil.o sgcond.o sgcharmx.o sggrid.o sgasym.o
	rm -f libsginfo.a
	ar r libsginfo.a sgvers.o sgclib.o sgio.o sgfind.o sghkl.o sgss.o sgharker.o sgutil.o sgcond.o sgcharmx.o sggrid.o sgasym.o

.c.o:
	$(CC) $(CCFLAGS) -c $*.c

all: sginfo sgtest sgquick libsginfo sgdev

clean:
	rm -f sgvers.o sgclib.o
	rm -f sgio.o sgfind.o sghkl.o sgss.o sgharker.o sgutil.o
	rm -f sginfo sginfo.o sgssvfy.o xtal32ss.o
	rm -f sgtest sgtest.o
	rm -f sgquick sgquick.o
	rm -f libsginfo.a core
	rm -f sgcond.o sgasciio.o sghklvfy.o
	rm -f sgdev.o sgdev sgcharmx.o sggrid.o sgasym.o cnslib.o sgnorm.o
