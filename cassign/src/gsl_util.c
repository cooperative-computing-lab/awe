/**
 * gsl_util.h
 *
 * Utilities for managing and accessing gsl_matrix types.
 */

#include "gsl_util.h"


/**
 * Transform a gsl_matrix into a single gsl_vector.
 * 
 * Parameters:
 *     mat - a populated matrix to flatten
 *     vec - the structure in which the flattened matrix will be stored
 *
 * Returns:
 *     An exit code representing the success of the flattening operation
 */
exit_t gsl_matrix_flatten (const gsl_matrix* mat, gsl_vector** vec) {
  // Allocate a vector that can contein every element of mat
  size_t N = mat->size1 * mat->size2;
  *vec = gsl_vector_calloc(N);
  
  int i, j, k=0;

  // Iterate over every row of the matrix
  for (i=0; i<mat->size1; i++) {

    // Iterate over each element of matrix row i
    for (j=0; j<mat->size2; j++) {
      
      // Get the value and store it at vec index k
      const double val = gsl_matrix_get(mat, i, j);
      gsl_vector_set(*vec, k, val);

      // k persists across row and column iterations
      k++;
    }
  }

  return exitOK;
}


/**
 * Print the contents of a gsl_vector to stdout.
 *
 * Parameters:
 *     vec - the vector to print
 *
 * Returns:
 *     void
 */
void gsl_vector_printf (const gsl_vector* vec) {
  printf ("[ ");
  
  // Iterate over each vector element and print it
  for (int i=0; i<vec->size; i++) {
    const double val = gsl_vector_get(vec, i);
    printf("%.3f ", val);
  }

  printf ("]");
}


/**
 * Print the conents of a gsl_matrix to stdout.
 *
 * Parameters:
 *     mat - the matrix to print
 *
 * Returns:
 *     void
 */
void gsl_matrix_printf (const gsl_matrix* mat) {
  printf ("{\n");

  // Group each row by brackets
  for (int r=0; r<mat->size1; r++) {
    printf ("%3s[", " ");

    // Print each element in the row
    for (int c=0; c<mat->size2; c++) {
      printf ("%.3f ", gsl_matrix_get (mat, r, c));
    }

    printf ("]\n");
  }

  printf ("}\n");
}


/**
 * Get a set of rows from a matrix and store them in a new matrix.
 *
 * Parameters:
 *     mat     - the matrix from which to get rows
 *     indices - the rows to get from the matrix
 *     result  - the matrix in which the rows will be stored
 *
 * Returns:
 *     An exit code representing the success of getting the rows
 */
exit_t gsl_matrix_get_rows (const gsl_matrix *mat, const gsl_vector *indices,
                                                        gsl_matrix **result) {
  // Ensure that the number of indices to get does not exceed the matrix rows
  assert (indices->size <= mat->size1);

  // Allocate a matrix of the appropriate size (see GSL documentation)
  *result = gsl_matrix_alloc (indices->size, mat->size2);

  // Iterate over the indices vector and get each matrix row index
  for (int i=0; i<indices->size; i++) {
    const size_t r = gsl_vector_get (indices, i);

    // Get each value in the specified matrix row and store it in result
    for (int c=0; c<mat->size2; c++) {
      const double val = gsl_matrix_get (mat, r, c);
      gsl_matrix_set (*result, i, c, val);
    }
  }

  return exitOK;
}
