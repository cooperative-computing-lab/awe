
#include "rmsd_calc.h"



theodata* theodata_alloc () {
  return (theodata*) malloc (sizeof(theodata));
}


exit_t theodata_init (const size_t nreal, const size_t ndim, theodata** theo) {
  const size_t npad = 4 + nreal - nreal % 4;

  *theo = theodata_alloc ();
  (*theo)->nreal = nreal;
  (*theo)->npad  = npad;
  (*theo)->ndim  = ndim;
  (*theo)->coords = gsl_matrix_calloc (ndim, npad);

  return exitOK;
}

void theodata_printf (const theodata* theo) {
  printf ("<theodata: nreal = %lu npad = %lu ndim = %lu G = %.3f>",
	  theo->nreal, theo->npad, theo->ndim, theo->g);
}


exit_t prepare_data (const gsl_matrix* mat, theodata** theo) {
  gsl_matrix* mat2 = gsl_matrix_calloc (mat->size1, mat->size2);
  gsl_matrix_memcpy (mat2, mat);
  center_structure (mat2);

  theodata_init (mat2->size1, mat2->size2, theo);
  center_structure (mat2);
  const double G = calculate_theo_g (mat2);
  (*theo)->g = G;

  for (int r=0; r<mat2->size1; r++) {
    for (int c=0; c<mat2->size2; c++) {
      const double v = gsl_matrix_get (mat2, r, c);
      gsl_matrix_set ((*theo)->coords, c, r, v);
    }}

  gsl_matrix_free (mat2);

  return exitOK;
}

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

  
