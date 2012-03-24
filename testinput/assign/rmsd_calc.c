
#include "rmsd_calc.h"



exit_t center_structure (gsl_matrix* mat) {

  double mean;
  double oldval, newval;
  gsl_vector_view view;
  for (int i=0; i<mat->size1; i++) {
    view = gsl_matrix_row (mat, i);
    if (view.vector.data == NULL) { return exitIX_OUT_OF_BOUNDS; }
    mean = gsl_stats_mean (view.vector.data, view.vector.stride, view.vector.size);
    for (int j=0; j<mat->size2; j++) {
      oldval = gsl_matrix_get (mat, i, j);
      newval = oldval - mean;
      gsl_matrix_set (mat, i, j, newval);
    }
  }

  return exitOK;
}

