
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

  celldata_printf (cell_data);

  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  xdrframe_printf (frame);

  printf ("~> Computing rmsds...\n");
  int assignment = 0;
  double minrmsd = DBL_MAX;
  int asns = 0;
  for (int c=0; c<cell_data->ncells; c++) {
    const gsl_matrix cell = celldata_get_cell (cell_data, c);
    const double rmsd = compute_rmsd (&cell, frame->coords);
    printf ("~~> rmsd[%5d]: %8.5f\n", c, rmsd);
    if (rmsd < minrmsd) {
      minrmsd = rmsd;
      assignment = c;
      asns ++;
    } }
  printf ("~> Assignment: %d jumps: %d\n", assignment, asns);


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
