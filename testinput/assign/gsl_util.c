
#include "gsl_util.h"

exit_t gsl_matrix_flatten (const gsl_matrix* mat, gsl_vector** vec) {
  size_t N = mat->size1 * mat->size2;
  *vec = gsl_vector_calloc (N);
  int k=0;

  for (int i=0; i<mat->size1; i++) {
    for (int j=0; j<mat->size2; j++) {
      const double val = gsl_matrix_get (mat, i, j);
      gsl_vector_set (*vec, k, val);
      k++;
    }}

  return exitOK;
}

void gsl_vector_printf (const gsl_vector* vec) {
  printf ("[ ");
  for (int i=0; i<vec->size; i++) {
    const double val = gsl_vector_get (vec, i);
    printf ("%.3f ", val);
  }
  printf ("]");
}

void gsl_matrix_printf (const gsl_matrix* mat) {
  printf ("{\n");
  for (int r=0; r<mat->size1; r++) {
    printf ("%3s[", " ");
    for (int c=0; c<mat->size2; c++) {
      printf ("%.3f ", gsl_matrix_get (mat, r, c));
    }
    printf ("]\n");
  }
  printf ("}\n");
}
