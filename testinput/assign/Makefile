

CC = clang
CPU = core2
CFLAGS = --std=c99 -I/opt/macports/current/include/ -g
# CFLAGS = --std=c99 -O2 -mtune=$(CPU) -fPIC -msse2 -msse3 -fopenmp -I/opt/macports/current/include/ -g

LDFLAGS = -L/opt/macports/current/lib -lxdrfile -lgsl -lm -framework veclib

CODE = assign.c celldata.c xdr_util.c exit_codes.c rmsd_calc.c gsl_util.c theobald_rmsd.c
HEADERS = $(CODE:.c:.h)
OBJS = $(CODE:.c=.o)

COMPILE = $(CC) $(CFLAGS) -c


a.out : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

assign.o : assign.c assign.h
	$(COMPILE) assign.c

celldata.o : celldata.c celldata.h
	$(COMPILE) celldata.c

xdr_util.o : xdr_util.c xdr_util.h
	$(COMPILE) xdr_util.c

exit_codes.o : exit_codes.c exit_codes.h
	$(COMPILE) exit_codes.c

rmsd_calc.o : rmsd_calc.c rmsd_calc.h
	$(COMPILE) rmsd_calc.c

theobald_rmsd.o : theobald_rmsd.c theobald_rmsd.h
	$(COMPILE) theobald_rmsd.c

gsl_util.o : gsl_util.c gsl_util.h
	$(COMPILE) gsl_util.c

clean :
	rm *.o a.out