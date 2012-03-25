
#include "assign.h"


int main (void) {

  const char* cells_file = "Gens.dat";
  const char* xtc_file   = "traj.xtc";

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
    /* printf ("~> Cell %d: ", c); */
    /* gsl_matrix_printf (&cell); */
    const double rmsd = compute_rmsd (&cell, frame->coords);
    printf ("~> rmsd[%d]: %.3f\n", c, rmsd);
    if (rmsd < minrmsd) {
      minrmsd = rmsd;
      assignment = c;
    }
  }

  printf ("~> Assignment: %d\n", assignment);



  exit (EXIT_SUCCESS);
}
