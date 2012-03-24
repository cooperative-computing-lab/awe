

CC = clang
CPU = core2
CFLAGS = --std=c99 -O2 -mtune=$(CPU) -fPIC -msse2 -msse3 -fopenmp -I/opt/macports/current/include/

LDFLAGS = -L/opt/macports/current/lib -lxdrfile -lgsl -lm

CODE = assign.c theobald_rmsd.c
OBJS = $(CODE:.c=.o)


a.out : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

$(OBJS) : $(CODE)
	$(CC) $(CODE) -c $(CFLAGS)

clean :
	rm *.o a.out