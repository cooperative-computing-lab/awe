
#include "assign.h"


int main (int argc, char *argv[]) {

  const char
    *cells_file = argv[1], // "Gens.dat";
    *xtc_file   = argv[2], // "traj.xtc";
    *out_file   = argv[3]; // assignment.dat

  printf ("~> Cells file: %s\n", cells_file);
  printf ("~> Xtc file: %s\n", xtc_file);

  celldata* cell_data;
  xdrframe* frame;

  celldata_load_file (cells_file, &cell_data);
  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  xdrframe_printf (frame);

  int assignment;
  double minrmsd = DBL_MAX;
  for (int c=0; c<cell_data->ncells; c++) {
    const gsl_matrix cell = celldata_get_cell (cell_data, c);
    const double rmsd = naive_3d_rmsd (&cell, frame->coords);
    printf ("~> rmsd[%d]: %.3f\n", c, rmsd);
    if (rmsd < minrmsd) {
      minrmsd = rmsd;
      assignment = c;
    }
  }

  printf ("~> Assignment: %d\n", assignment);

  FILE *out = fopen (out_file, "w");
  if (out == NULL) {
    char emsg[50];
    sprintf (emsg, "Could not open for writing: %s\n", out_file);
    perror (emsg);
    exit (EXIT_FAILURE);
  }

  fprintf (out, "%d", assignment);
  fclose (out);

  exit (EXIT_SUCCESS);
}
