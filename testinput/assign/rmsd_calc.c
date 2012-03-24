
#include "rmsd_calc.h"


// TODO: center_structure should 
exit_t center_structure (gsl_matrix* mat) {

  // compute the mean along each axis (column)
  gsl_vector means = *gsl_vector_calloc (mat->size2);
  for (int d=0; d<mat->size2; d++) {
    const gsl_vector_const_view view = gsl_matrix_const_column (mat, d);
    const double mean = gsl_stats_mean (view.vector.data, view.vector.stride, view.vector.size);
    gsl_vector_set (&means, d, mean);
  }

  // center the coordinates
  for (int r=0; r<mat->size1; r++) {
    for (int d=0; d<mat->size2; d++) {
      const double m = gsl_vector_get (&means, d);
      const double v0 = gsl_matrix_get (mat, r, d);
      const double v1 = v0 - m;
      gsl_matrix_set (mat, r, d, v1);
    }}

  return exitOK;
}


double calculate_theo_g (const gsl_matrix* mat) {

  double ssqrs = 0;
  for (int d=0; d<mat->size2; d++) {
    const gsl_vector_const_view col = gsl_matrix_const_column (mat, d);
    for (int i=0; i<col.vector.size; i++) {
      double s = gsl_vector_get (&col.vector, i);
      s *= s;

      // sum of squares
      ssqrs += s;
    }}

  return ssqrs;
}
