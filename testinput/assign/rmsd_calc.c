
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
  (*theo)->coords = (float*) calloc (ndim * npad, sizeof(float));

  return exitOK;
}

void theodata_printf (const theodata* theo) {
  printf ("<theodata: nreal = %lu npad = %lu ndim = %lu G = %.3f>",
	  theo->nreal, theo->npad, theo->ndim, theo->g);
}

void theodata_printf_all (const theodata* theo) {
  theodata_printf (theo);
  printf ("{\n");
  for (int r=0; r<theo->ndim; r++) {
    printf ("%3s[", " ");
    for (int c=0; c<theo->npad; c++) {
      printf ("%.3f ", theodata_get (theo, r, c));
    }
    printf ("]\n");
  }
  printf ("}\n");
}

exit_t prepare_data (const gsl_matrix* mat, theodata** theo) {
  size_t
    Rs = mat->size1,
    Cs = mat->size2;

  gsl_matrix* mat2 = gsl_matrix_calloc (Rs, Cs);
  theodata_init (Rs, Cs, theo);


  gsl_matrix_memcpy (mat2, mat);
  center_structure (mat2);
  /* printf ("~> Centered coordinates: "); */
  /* gsl_matrix_printf (mat2); */

  const double G = calculate_theo_g (mat2);
  (*theo)->g = G;

  for (int r=0; r<Rs; r++) {
    for (int c=0; c<Cs; c++) {
      const float v = (float) gsl_matrix_get (mat2, r, c);
      theodata_set ((*theo), c, r, v);
    }}

  gsl_matrix_free (mat2);

  return exitOK;
}

size_t theodata_ix (const theodata* theo, const size_t dim, const size_t coord) {
  return dim * theo->npad + coord;
}

exit_t theodata_set (theodata* theo, const size_t dim, const size_t coord, const float val) {
  const size_t ix = theodata_ix (theo, dim, coord);
  theo->coords[ix] = val;
  return exitOK;
}

float theodata_get (const theodata* theo, const size_t dim, const size_t coord) {
  const size_t ix = theodata_ix (theo, dim, coord);
  return theo->coords[ix];
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

double theo_rmsd (const theodata* theo1, const theodata* theo2) {

  assert (theo1->nreal	== theo2->nreal);
  assert (theo1->npad	== theo2->npad);
  assert (theo1->ndim	== theo2->ndim);

  return ls_rmsd2_aligned_T_g (theo1->nreal, theo1->npad, theo1->npad, theo1->coords, theo2->coords, theo1->g, theo2->g);
}


double compute_rmsd (const gsl_matrix* m1, const gsl_matrix* m2) {
  theodata *theo1, *theo2;

  prepare_data (m1, &theo1);
  prepare_data (m2, &theo2);
  return theo_rmsd (theo1, theo2);
}
