# Unix-Makefile for kriber

CC=cc
CCFLAGS = -O
CCWARN  =

kriber.x: kriber.o p2clib.o
	$(CC) $(CCFLAGS) -o kriber.x kriber.o p2clib.o -lm
#	strip kriber.x

kriber.o: kriber.c p2c.h
	$(CC) $(CCFLAGS) $(CCWARN) -c kriber.c

p2clib.o: p2clib.c p2c.h
	$(CC) $(CCFLAGS) -c p2clib.c

clean:
	rm -f kriber.o p2clib.o kriber.x
