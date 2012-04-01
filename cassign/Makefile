

CC = gcc
CFLAGS = --std=c99 -g
LDFLAGS = -lxdrfile -lgsl -lgslcblas -lm

CODE = assign.c celldata.c xdr_util.c exit_codes.c rmsd_calc.c gsl_util.c theobald_rmsd.c
HEADERS = $(CODE:.c:.h)
OBJS = $(CODE:.c=.o)

COMPILE = $(CC) $(CFLAGS) -c



assign : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o assign

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
	rm *.o assign

