

CC = clang
CPU = core2
CFLAGS = --std=c99 -I/opt/macports/current/include/ -g
# CFLAGS = --std=c99 -O2 -mtune=$(CPU) -fPIC -msse2 -msse3 -fopenmp -I/opt/macports/current/include/ -g

LDFLAGS = -L/opt/macports/current/lib -lxdrfile -lgsl -lm

CODE = assign.c celldata.c xdr_util.c theobald_rmsd.c
HEADERS = $(CODE:.c:.h)
OBJS = $(CODE:.c=.o)


a.out : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

assign.o : assign.c assign.h
	$(CC) -c assign.c $(CFLAGS)

theobald_rmsd.o : theobald_rmsd.c theobald_rmsd.h
	$(CC) -c theobald_rmsd.c $(CFLAGS)

clean :
	rm *.o a.out