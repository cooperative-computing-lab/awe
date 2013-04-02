

CC = gcc
CFLAGS = --std=c99 -g
LDFLAGS = -lxdrfile -lgsl -lgslcblas -lm # -static

CODE = $(addprefix src/, assign.c celldata.c xdr_util.c exit_codes.c rmsd_calc.c gsl_util.c)
HEADERS = $(CODE:.c:.h)
OBJS = $(CODE:.c=.o)

COMPILE = $(CC) $(CFLAGS) -c


TESTFILES = $(addprefix testfiles/, Gens.dat traj.xtc AtomIndices.dat)



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
	rm -f src/*.o assign

check : assign $(TESTFILES)
	$< testfiles/Gens.dat testfiles/AtomIndices.dat testfiles/traj.xtc testfiles/AtomIndices.dat testfiles/cellasn.dat
