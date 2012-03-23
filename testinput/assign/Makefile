

CC= clang
CPU = core2
CFLAGS = --std=c99 -O2 -mtune=$(CPU) -fPIC -shared -msse2 -msse3 -fopenmp

LDFLAGS = -lxdrfile

a.out : load-datfile.o
	$(CC) load-datfile.o $(LDFLAGS)

load-datfile.o : load-datfile.c
	$(CC) -c load-datfile.c

clean :
	rm -v *.o a.out