

CC= clang
CPU = core2
CFLAGS = --std=c99 -O2 -mtune=$(CPU) -fPIC -shared -msse2 -msse3 -fopenmp

LDFLAGS = -lxdrfile

a.out : assign.o
	$(CC) assign.o $(LDFLAGS)

assign.o : assign.c
	$(CC) -c assign.c

clean :
	rm -v *.o a.out