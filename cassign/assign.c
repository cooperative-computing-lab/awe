
#include "assign.h"


int main (int argc, char *argv[]) {

  const char
    *cells_file = /*argv[1], */ "Gens.dat",
    *xtc_file   = /*argv[2], */ "traj.xtc",
    *mndx_file  = /*argv[3], */ "AtomIndices.dat",
    *out_file   = /*argv[4]; */ "asn.dat";

  printf ("~> Cells file: %s\n", cells_file);
  printf ("~> Xtc file: %s\n", xtc_file);

  celldata* cell_data;
  xdrframe* frame;
  gsl_vector *atomindices;

  celldata_load_file (cells_file, &cell_data);

  celldata_printf (cell_data);

  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  xdrframe_printf (frame);

  xdrframe_load_atomindices (mndx_file, &atomindices);
  gsl_vector_printf (atomindices);

  celldata *newcells;
  xdrframe *newframe;
  celldata_get_rows (cell_data, atomindices, &newcells);
  xdrframe_select_atoms (frame, atomindices, &newframe);

  celldata_printf (newcells);
  xdrframe_printf (newframe);

  

  printf ("~> Computing rmsds...\n");
  int assignment = 0;
  double
    minrmsd = DBL_MAX,
    maxrmsd = -1;
  int asns = 0;
  for (int c=0; c<cell_data->ncells; c++) {
    const gsl_matrix cell = celldata_get_cell (newcells, c);
    const double rmsd = compute_rmsd (&cell, newframe->coords);
    printf ("~~> rmsd[%5d]: %8.5f\n", c, rmsd);
    if (rmsd < minrmsd) {
      minrmsd = rmsd;
      assignment = c;
      asns ++;
    }
    if (rmsd > maxrmsd ) { maxrmsd = rmsd; }
  }
  printf ("~> Assignment: %d jumps: %d\n", assignment, asns);
  printf ("~> minrmsd: %f maxrmsd: %f\n", minrmsd, maxrmsd);


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
