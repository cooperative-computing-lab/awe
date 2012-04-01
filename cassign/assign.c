
#include "assign.h"

exit_t load_atomindices (const char *mndxpath, gsl_vector **target) {
  /** 1) get the number of lines the file
      2) allocate the vector
      3) reread the file, inserting values into the vector

      This access the file twice, buf unfortunately, I'm not (yet)
      aware of any standard list datastructure for C, and am willing
      to accept the redundancy at the moment
  */

  exit_t status = exitOK;

  int linecount = 0;
  size_t ndx;

  FILE *mndx = fopen (mndxpath, "r");
  if (mndx == NULL) { return exitPATH_NOT_FOUND; }

  while ( ! feof (mndx) ) {
    int count = fscanf (mndx, "%d", &ndx);
    if (count == 1)
      { linecount ++; }
  }

  fclose (mndx);
  if ( ! status == exitOK ) { return status; }

  mndx = fopen (mndxpath, "r");

  *target = gsl_vector_calloc (linecount);
  size_t i = 0;

  while ( ! feof (mndx) ) {
    int count = fscanf (mndx, "%d", &ndx);
    if (count == 1) {
      gsl_vector_set (*target, i, (int)ndx);
      i++;
    }
    else if ( count > 1 ) {
      perror ("Atom indices file incorrectly formatted");
      status = exitFAILURE;
      break;
    }
    else { break; }
  }

  fclose (mndx);

  return status;
}



int main (int argc, char *argv[]) {

  const char
    *cells_file	   = argv[1], //"Gens.dat",
    *cell_ndx_file = argv[2], //"AtomIndices.dat",
    *xtc_file	   = argv[3], //"traj.xtc",
    *xtc_ndx_file  = argv[4], //"AtomIndices.dat",
    *out_file	   = argv[5]; //"cell2.dat";

  printf ("~> Cells file: %s\n", cells_file);
  printf ("~> Xtc file: %s\n", xtc_file);

  celldata* cell_data;
  xdrframe* frame;
  gsl_vector *cell_indices, *xtc_indices;

  celldata_load_file (cells_file, &cell_data);
  load_atomindices (cell_ndx_file, &cell_indices);

  printf ("~> Loaded cells: ");
  celldata_printinfo (cell_data);
  printf ("\n");

  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  printf ("\n");
  /* xdrframe_printf (frame); */

  load_atomindices (xtc_ndx_file, &xtc_indices);
  printf ("~> Atom Indices: ");
  gsl_vector_printf (xtc_indices);
  printf ("\n");

  celldata *newcells;
  xdrframe *newframe;
  celldata_get_rows (cell_data, cell_indices, &newcells);
  xdrframe_select_atoms (frame, xtc_indices , &newframe);

  if (newcells->ncoords != newframe->natoms) {
    printf ("Number of atoms mismatch after selecting atom indices: cell = %lu structure = %lu\n", newcells->ncoords, newframe->natoms);
    exit (1);
  }

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
