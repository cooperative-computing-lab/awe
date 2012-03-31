
#include "assign.h"


int main (int argc, char *argv[]) {

  const char
    *cells_file = argv[1], //"Gens.dat",
    *xtc_file   = argv[2], //"traj.xtc",
    *mndx_file  = argv[3], //"AtomIndices.dat",
    *out_file   = argv[4]; //"cell2.dat";

  printf ("~> Cells file: %s\n", cells_file);
  printf ("~> Xtc file: %s\n", xtc_file);

  celldata* cell_data;
  xdrframe* frame;
  gsl_vector *atomindices;

  celldata_load_file (cells_file, &cell_data);

  printf ("~> Loaded cells: ");
  celldata_printinfo (cell_data);
  printf ("\n");

  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  printf ("\n");
  /* xdrframe_printf (frame); */

  xdrframe_load_atomindices (mndx_file, &atomindices);
  printf ("~> Atom Indices: ");
  gsl_vector_printf (atomindices);
  printf ("\n");

  celldata *newcells;
  xdrframe *newframe;
  celldata_get_rows (cell_data, atomindices, &newcells);
  xdrframe_select_atoms (frame, atomindices, &newframe);

  printf ("~> Using: ");
  celldata_printinfo (newcells);
  printf ("\n");
  printf ("~> Using: ");
  xdrframe_printsummary (newframe);
  printf ("\n");
  

  printf ("~> Computing rmsds...\n");
  int assignment = 0;
  double
    minrmsd = DBL_MAX,
    maxrmsd = -1;
  int asns = 0;
  for (int c=0; c<cell_data->ncells; c++) {
    const gsl_matrix cell = celldata_get_cell (newcells, c);
    const gsl_matrix *coords = newframe->coords;
    const double rmsd = compute_rmsd (&cell, coords);
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


  printf ("~> Saving assignment to: %s\n", out_file);
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
