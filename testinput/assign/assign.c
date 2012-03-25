
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

  theodata* theo;
  prepare_data (frame->coords, &theo);
  printf ("~> Prepared data for frames:\n");
  theodata_printf_all (theo);


  exit (EXIT_SUCCESS);
}
