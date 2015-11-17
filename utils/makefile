CC=cc
CCFLAGS = -O
CCWARN  =

all:
	$(CC) $(CCFLAGS) section.c -o section
	$(CC) $(CCFLAGS) genseq.c -o genseq
	$(CC) $(CCFLAGS) coseq_reduce.c -o coseq_reduce
	$(CC) $(CCFLAGS) coseq_cmp.c -o coseq_cmp
	$(CC) $(CCFLAGS) xrs2scandata.c -o xrs2scandata

install:
	mv coseq_cmp ~/bin/zeolites/
	mv coseq_reduce ~/bin/zeolites/
	mv genseq ~/bin/zeolites/
	mv xrs2scandata ~/bin/zeolites/
	mv section ~/bin/zeolites/

clean:
	rm -f coseq_cmp
	rm -f coseq_reduce
	rm -f genseq
	rm -f xrs2scandata
	rm -f section