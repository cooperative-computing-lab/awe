/**
 * rmsd_calc.c
 *
 * Functions for calculating the Root Mean Square Distance between two
 * sets of coordinates (e.g. molecular atomic coordinates).
 */

#include "rmsd_calc.h"


/**
 * Normalize coordinates using the mean dimension values.
 *
 * Parameters:
 *     mat - the matrix to center
 *
 * Returns:
 *     
 */
exit_t center_structure(gsl_matrix* mat) {
  // There seems to be a misunderstanding of what 'const' means and a lot of
  // unnecessary variable initialization.

  // compute the mean along each axis (column)
  gsl_vector means = *gsl_vector_calloc(mat->size2);
  for (int d=0; d<mat->size2; d++) {
    const gsl_vector_const_view view = gsl_matrix_const_column(mat, d);
    const double mean = gsl_stats_mean(view.vector.data, view.vector.stride, view.vector.size);
    gsl_vector_set(&means, d, mean);
  }

  // Fix the unnecessary initializations
  // center the coordinates
  for (int r=0; r<mat->size1; r++) {
    for (int d=0; d<mat->size2; d++) {
      const double m = gsl_vector_get (&means, d);
      const double v0 = gsl_matrix_get (mat, r, d);
      const double v1 = v0 - m;
      gsl_matrix_set (mat, r, d, v1);
    }
  }

  return exitOK;
}


/**
 * Alias for the function computing the RMSD for ease of refactoring.
 * The result should be the same in either order.
 *
 * Parameters:
 *     m1 - a reference matrix
 *     m2 - a reference matrix of the same shape as m1
 *
 * Returns:
 *     A double representing the RMSD between m1 and m2
 */
double compute_rmsd (const gsl_matrix* m1, const gsl_matrix* m2) {
  return kabsch_rmsd (m1, m2);
}


/**
 * Prepare the matrices and them compute the RMSD.
 *
 * Parameters:
 *     m1 - a reference matrix
 *     m2 - a reference matrix of the same shape as m1
 *
 * Returns:
 *     A double representing the RMSD between m1 and m2
 */
double kabsch_rmsd (const gsl_matrix *m1, const gsl_matrix *m2) {
  // Ensure that the matrices have the same shape
  assert (m1->size1 == m2->size1);
  assert (m1->size2 == m2->size2);

  // Get the shape of the matrices
  const int
    N = m1->size1,
    D = m1->size2;
    

  // Get the transposes of m1 and m2
  gsl_matrix
    *P = gsl_matrix_calloc (D, N),
    *Q = gsl_matrix_calloc (D, N);

  gsl_matrix_transpose_memcpy (P, m1);
  gsl_matrix_transpose_memcpy (Q, m2);


  // Containers for returns from the RMSD function
  gsl_matrix *U, *t;
  double rmsd;
  kabsch_function (P, Q, &U, &t, &rmsd);

  // Free up the allocated memory
  gsl_matrix_free (P);
  gsl_matrix_free (Q);
  gsl_matrix_free (U);
  gsl_matrix_free (t);

  return rmsd;

}  


/**
 * Find the Least Root Mean Square between two sets of N points in D dimensions
 * and the rigid transformation (i.e. translation and rotation) 
 * to employ in order to bring one set that close to the other,
 * using the Kabsch (1976) algorithm.
 *
 * Parameters
 *     P - a D*N matrix where P(a,i) is the a-th coordinate of the i-th point 
 *         in the 1st representation
 *     Q - a D*N matrix where Q(a,i) is the a-th coordinate of the i-th point 
 *         in the 2nd representation
 *
 * Output
 *     U    - a proper orthogonal D*D matrix, representing the rotation
 *     r    - a D-dimensional column vector, representing the translation
 *     rmsd - Root Mean Square Distance
 *
 * References:
 *   1) Kabsch W. A solution for the best rotation to relate two sets of vectors. Acta Cryst A 1976;32:9223.
 *   2) Kabsch W. A discussion of the solution for the best rotation to relate two sets of vectors. Acta Cryst A 1978;34:8278.
 *   3) http://cnx.org/content/m11608/latest/
 *   4) http://en.wikipedia.org/wiki/Kabsch_algorithm
 */
exit_t kabsch_function (const gsl_matrix *A, const gsl_matrix *B, gsl_matrix **U, gsl_matrix **r, double *rmsd) {

  assert (A->size1 == B->size1);
  assert (A->size2 == B->size2);

  const int
    D = B->size1,
    N = A->size2;


  gsl_matrix
    *P = gsl_matrix_calloc (D, N),
    *Q = gsl_matrix_calloc (D, N);

  gsl_matrix_memcpy (P, A);
  gsl_matrix_memcpy (Q, B);



  // m = ones(1,N)/N
  gsl_matrix *m = gsl_matrix_calloc (1, N);
  gsl_matrix_set_all (m, 1);
  gsl_matrix_scale (m, 1.0/(double)N);
  /* OK: [ 0.045 0.045 ... ]
   */


  // p0 = P * m' // centroid of P
  gsl_matrix *p0 = gsl_matrix_calloc (D,1);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, P, m, 0.0, p0);
  /* OK: [2.304 ]
         [0.896 ]
         [-0.181 ]
  */


  // q0 = Q * m' // centroid of Q
  gsl_matrix *q0 = gsl_matrix_calloc (D, 1);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Q, m, 0.0, q0);
  /* OK: [-5.174 ]
         [5.144 ]
         [9.810 ]
  */


  // v1 = ones(1,N)
  gsl_matrix *v1 = gsl_matrix_calloc (1, N);
  gsl_matrix_set_all (v1, 1.0);
  /* OK: [1.0000 ... ]
   */


  // P = P - p0*v1 ;      % translating P to center the origin
  gsl_matrix *mtmp = gsl_matrix_calloc (D, N);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, p0, v1, 0.0, mtmp);
  /* OK:    [2.304 ... ]
            [0.896 ... ]
            [-0.181 ...]
  */
  gsl_matrix_sub (P, mtmp);
  /* OK :   [-0.259 -0.164 -0.153 -0.152 -0.050 0.036 -0.055 -0.129 0.043 0.059 0.162 0.199 0.128 0.248 -0.013 -0.128 0.060 0.144 0.011 -0.068 0.088 -0.014 ]
            [-0.190 -0.191 -0.117 -0.286 -0.167 -0.252 -0.050 0.022 -0.021 -0.117 0.055 -0.013 0.147 0.069 0.075 0.125 0.099 0.051 0.176 0.128 0.183 0.277 ]
        [-0.222 -0.274 -0.357 -0.322 -0.175 -0.161 -0.111 -0.143 -0.005 0.037 -0.066 -0.147 -0.107 0.002 0.102 0.091 0.217 0.238 0.330 0.385 0.401 0.293 ]
  */
  gsl_matrix_free (mtmp);


  // Q = Q - q0*v1 ;      % translating Q to center the origin
  mtmp = gsl_matrix_calloc (D, N);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, q0, v1, 0.0, mtmp);
  /* OK:    [2.304 ... ]
            [0.896 ... ]
            [-0.181 ...]
  */
  gsl_matrix_sub (Q, mtmp);
  /* OK :   [0.405 0.326 0.324 0.332 0.201 0.201 0.097 0.109 -0.028 -0.005 -0.081 -0.004 -0.118 -0.164 -0.125 -0.163 -0.155 -0.130 -0.238 -0.243 -0.341 -0.205 ]
            [-0.008 0.046 0.015 0.154 0.003 -0.092 0.084 0.172 0.074 0.044 0.216 0.289 0.237 0.212 -0.027 -0.017 -0.144 -0.138 -0.258 -0.322 -0.231 -0.309 ]
            [-0.045 -0.095 -0.200 -0.089 -0.018 0.058 -0.028 -0.076 0.044 0.146 0.055 0.079 -0.046 0.126 -0.020 -0.138 0.049 0.146 0.021 0.109 -0.000 -0.070 ]
  */
  gsl_matrix_free (mtmp);




  // C is a covariance matrix of the coordinates
  // C = P*Q' / N
  gsl_matrix *C = gsl_matrix_calloc (D, D);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, P, Q, 0.0, C);
  gsl_matrix_scale (C, 1.0 / (double) N);
  /* OK:    [-0.017 0.004 0.008 ]
            [-0.029 -0.010 0.002 ]
            [-0.043 -0.028 0.009 ]
  */


  // [V,S,W] = svd(C)
  gsl_matrix *W     = gsl_matrix_calloc (D, D);
  gsl_vector *s     = gsl_vector_calloc (D);
  gsl_vector *work  = gsl_vector_calloc (D);
  gsl_matrix *V     = gsl_matrix_calloc (D, D);       // SV_decomp would overwrite C to store U
  gsl_matrix_memcpy             (V, C);
  gsl_linalg_SV_decomp    (V, W, s, work);
  /* OK:
     V:    [-0.230 -0.903 -0.364 ]
           [-0.489 -0.216 0.845 ]
           [-0.841 0.372 -0.392 ]

     s:    [ 0.062 0.014 0.005 ]

     W:    [0.876 0.387 -0.289 ]
           [0.453 -0.866 0.210 ]
           [-0.169 -0.315 -0.934 ]
  */




  // d = sign(det(C))
  gsl_matrix *LU = gsl_matrix_calloc (D, D);
  gsl_permutation *perm = gsl_permutation_calloc (D);
  int signum = 1;
  gsl_matrix_memcpy (LU, C);
  gsl_linalg_LU_decomp (LU, perm, &signum);
  const double determinant = gsl_linalg_LU_det (LU, signum);
  const int d = determinant > 0 ? 1 : -1;
  gsl_matrix_free (LU);
  gsl_permutation_free (perm);
  /* OK: determinant: 4.669769e-06
         d:           1
  */



  // I      = eye(D)
  // I(D, D)    = sign(det(C))
  gsl_matrix *I = gsl_matrix_calloc (D, D);
  gsl_matrix_set_identity (I);
  gsl_matrix_set (I, D-1, D-1, d);
  /* OK: I:    [1.000 0.000 0.000 ]
               [0.000 1.000 0.000 ]
               [0.000 0.000 1.000 ]
  */
  



  // U = W*I*V'
  *U = gsl_matrix_calloc (D, D);
  mtmp = gsl_matrix_calloc (D,D);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, I, V, 0.0, mtmp);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, W, mtmp, 1.0, *U);
  /* OK: U:    [-0.446 -0.756 -0.479 ]
               [0.601 0.144 -0.786 ]
               [0.663 -0.638 0.391 ]
  */
  gsl_matrix_free (mtmp);




  // r : a D-dimensional column vector, representing the translation
  // r = q0 - U*p0
  *r            = gsl_matrix_calloc (D, 1);
  mtmp          = gsl_matrix_calloc (D, 1);
  gsl_matrix_memcpy       (*r, q0)                 ;                  // gsl_matrix_sub clobbers 1st parameter
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, *U, p0, 0.0, mtmp);
  gsl_matrix_sub (*r, mtmp);
  /* OK: r:    [-3.557 ]
               [3.489 ]
           [8.925 ]
  */
  gsl_matrix_free (mtmp);



  // diff = U*P - Q
  gsl_matrix *diff = gsl_matrix_calloc (D, N);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, *U, P, 1.0, diff);
  gsl_matrix_sub(diff, Q);
  /* OK:    [-0.040 0.023 0.004 0.106 0.031 0.051 0.018 0.000 0.027 0.049 -0.001 -0.005 0.001 0.000 0.025 0.082 -0.051 -0.087 -0.058 -0.008 -0.029 -0.139 ]
            [0.000 0.043 0.157 -0.034 0.080 0.205 -0.036 -0.134 -0.047 -0.054 -0.058 -0.055 -0.054 -0.055 -0.049 -0.114 0.024 0.045 0.031 -0.002 -0.004 0.110 ]
        [-0.092 0.002 0.034 0.046 0.024 0.064 -0.019 -0.079 -0.003 -0.017 -0.008 0.004 -0.004 -0.004 0.004 0.009 0.013 0.010 0.003 -0.085 0.099 -0.001 ]
  */


  //lrms = 0 ;
  //for i=1:N
  //    lrms = lrms + m(i)*Diff(:,i)'*Diff(:,i) ;
  //end
  double lrms = 0.0;
  double dot = 0.0;
  for (int i=0; i<N; i++) {
    gsl_vector_view xyz = gsl_matrix_column (diff, i);
    gsl_blas_ddot (&xyz.vector, &xyz.vector, &dot);
    lrms = lrms + gsl_matrix_get (m, 0, i) * dot;
  }

  *rmsd = sqrt (lrms);



  // clean up
  gsl_matrix_free (P);
  gsl_matrix_free (Q);
  gsl_matrix_free (m);
  gsl_matrix_free (p0);
  gsl_matrix_free (q0);
  gsl_matrix_free (v1);
  gsl_matrix_free (C);
  gsl_matrix_free (W);
  gsl_vector_free (s);
  gsl_vector_free (work);
  gsl_matrix_free (V);
  gsl_matrix_free (I);
  gsl_matrix_free (diff);


  return exitOK;


}
