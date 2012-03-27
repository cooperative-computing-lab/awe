
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
      printf ("%.3f, ", theodata_get (theo, r, c));
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

  double G = 0;
  for (int d=0; d<mat->size2; d++) {
    const gsl_vector_const_view col = gsl_matrix_const_column (mat, d);
    double dot = 0;
    gsl_blas_ddot (&col.vector, &col.vector, &dot);
    G += dot;
  }

    /* for (int i=0; i<col.vector.size; i++) { */
    /*   double s = gsl_vector_get (&col.vector, i); */
    /*   s *= s; */

    /*   // sum of squares */
    /*   ssqrs += s; */
    /* }} */

  return G; // ssqrs;
}

double theo_rmsd (const theodata* theo1, const theodata* theo2) {

  assert (theo1->nreal	== theo2->nreal);
  assert (theo1->npad	== theo2->npad);
  assert (theo1->ndim	== theo2->ndim);

  const double msd = ls_rmsd2_aligned_T_g (theo1->nreal, theo1->npad, theo1->npad, theo1->coords, theo2->coords, theo1->g, theo2->g);
  return sqrt (msd);
}


double compute_rmsd (const gsl_matrix* m1, const gsl_matrix* m2) {
  theodata *theo1, *theo2;

  prepare_data (m1, &theo1);
  prepare_data (m2, &theo2);

  /* theodata_printf_all (theo1); */
  /* theodata_printf_all (theo2); */

  return theo_rmsd (theo1, theo2);
}


/** DO NOT USE */
double naive_3d_rmsd (const gsl_matrix* m1, const gsl_matrix* m2) {

  /** http://cnx.org/content/m11608/latest/ */

  // sanity check
  assert (m1->size2 == 3);
  assert (m2->size2 == 3);
  assert (m1->size1 == m2->size1);

  const size_t D = 3;
  const size_t N = m1->size1;

  // 0) center structures
  gsl_matrix
    *x = gsl_matrix_alloc (N, D),
    *y = gsl_matrix_alloc (N, D);
  gsl_matrix_memcpy (x, m1);
  gsl_matrix_memcpy (y, m2);
  center_structure (x);
  center_structure (y);


  // 1) build 3XN matrices X and Y as transposed x, y
  gsl_matrix
    *X = gsl_matrix_alloc (D, N),
    *Y = gsl_matrix_alloc (D, N);
  gsl_matrix_transpose_memcpy (X, x);
  gsl_matrix_transpose_memcpy (Y, y);


  // 2) Compute the covariance matrix C = XY'
  gsl_matrix *C = gsl_matrix_calloc (D, D);
  gsl_blas_dgemm ( CblasNoTrans, CblasTrans, 1.0, X, Y, 0.0, C);

  // 3) Compute the SVD (Singular Value Decomposition) of C = VSW'
  gsl_matrix
    *V = gsl_matrix_calloc (D, D),
    *W = gsl_matrix_calloc (D, D);
  gsl_vector
    *s    = gsl_vector_calloc (D),
    *work = gsl_vector_calloc (D);

  // gsl_linalg_SV_decomp replaced the first matrix with V
  gsl_matrix_memcpy (V, C);
  gsl_linalg_SV_decomp (V, W, s, work);

  // 4) Compute d = sign(det(C))
  gsl_matrix *LU = gsl_matrix_calloc (D, D);
  gsl_permutation *perm = gsl_permutation_calloc (D);
  int signum = 1;
  gsl_matrix_memcpy (LU, C);
  gsl_linalg_LU_decomp (LU, perm, &signum);
  const int determinant = gsl_linalg_LU_det (LU, signum);
  const int d = determinant > 0 ? 1 : -1;
  gsl_matrix_free (LU);
  gsl_permutation_free (perm);

  // 5) Compute the optimal rotation U as
  /** 
            | 1 0 0 |
      U = W | 0 1 0 | V'
            | 0 0 d |
  */
  gsl_matrix *I = gsl_matrix_calloc (D, D);
  gsl_matrix_set_identity (I);
  gsl_matrix_set (I, D-1, D-1, d);
  gsl_matrix_transpose (V);
  gsl_matrix 
    *u = gsl_matrix_calloc (D, D),
    *U = gsl_matrix_calloc (D, D);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, W, I, 0.0, u);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, u, V, 0.0, U);


  // 6) fit Y to X
  gsl_matrix *Ycp = gsl_matrix_alloc (Y->size1, Y->size2);
  gsl_matrix_memcpy (Ycp, Y);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, U, Ycp, 0.0, Y);

  // Compute RMSD
  /** The RMSD is the square root of the average of the squared distances between corresponding atoms of X and Y.
            ______________________________
           /       N
          /       ____
         /     1  \    |         | 2
 _      /     ___  \   | X - Y   |
  \    /       N   /   |  i    i |
   \  /           /___
    \/	          i = 1
  */

  
  double sum_squares = 0.0;
  for (int i=0; i <N; i++) {
    gsl_vector_view x = gsl_matrix_column (X, i);
    gsl_vector_view y = gsl_matrix_column (Y, i);
    gsl_vector_sub (&x.vector, &y.vector);
    double dot = 0.0;
    gsl_blas_ddot (&x.vector, &x.vector, &dot);
    sum_squares += dot;
  }
  const double rmsd = sqrt (sum_squares / N);

  // free resources
  gsl_matrix_free (x);
  gsl_matrix_free (y);
  gsl_matrix_free (X);
  gsl_matrix_free (Y);
  gsl_matrix_free (C);
  gsl_matrix_free (V);
  gsl_matrix_free (W);
  gsl_vector_free (s);
  gsl_vector_free (work);
  gsl_matrix_free (I);
  gsl_matrix_free (u);
  gsl_matrix_free (U);
  gsl_matrix_free (Ycp);


  // done
  return rmsd;


}
