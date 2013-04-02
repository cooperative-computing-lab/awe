
#include "gsl_util.h"
#include "new_rmsd.h"


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <stdio.h>


void load_data (const char *path, const int N, const int D, gsl_matrix **dest) {

  *dest = gsl_matrix_calloc (N, D);

  float x, y, z;

  int atom = 0;
  FILE *fd = fopen (path, "r");
  while ( ! feof (fd) ) {
    if ( fscanf (fd, "%f %f %f", &x, &y, &z) == 3) {
      gsl_matrix_set (*dest, atom, 0, (double) x);
      gsl_matrix_set (*dest, atom, 1, (double) y);
      gsl_matrix_set (*dest, atom, 2, (double) z);
      atom++;
    }
  }

  fclose (fd);

}


int main (const int argc, const char *argv[]) {

  gsl_matrix *cell, *frame;
  float x, y, z;
  FILE *fd;

  const int N = 22, D = 3;
  cell  = gsl_matrix_calloc (N, D);
  frame = gsl_matrix_calloc (N, D);

  const char 
    *cellfile  = "cell0.dat",
    *framefile = "frame.dat";


  load_data (cellfile,  N, D, &cell);
  load_data (framefile, N, D, &frame);

  gsl_matrix_printf (cell);
  gsl_matrix_printf (frame);

  const double r = rmsd (cell, frame);
  printf ("rmsd: %f\n", r);

  exit (0);
}
