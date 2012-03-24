
#include "assign.h"







/* float compute_rmsd (float** ref, float** structure, int natoms, int ndims) { */

/*   float *AData, *BData; */
/*   return 42; */

/*   /\* float msd = ls_rmsd2_aligned_T_g(nrealatoms,npaddedatoms,rowstride,(AData+i*truestride),BData,GAData[i],G_y); *\/ */

/* } */
  


int main (void) {

  const char* cells_file = "Gens.dat";
  const char* xtc_file   = "traj.xtc";

  celldata* cell_data;
  xdrframe* frame;

  celldata_load_file (cells_file, cell_data);
  xdrframe_last_in_xtc (xtc_file, &frame);

  xdrframe_printsummary (frame);
  xdrframe_printf (frame);

  theodata* theo;
  prepare_data (frame->coords, &theo);
  printf ("~> Prepared data for frame:\n");
  theodata_printf (theo);
  gsl_matrix_printf (theo->coords);


  exit (EXIT_SUCCESS);
}
