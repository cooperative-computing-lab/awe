/**
 * gsl_util.h
 *
 * Utilities for managing and accessing gsl_matrix types.
 */

#ifndef _GSL_UTIL_H_
#define _GSL_UTIL_H_


#include "exit_codes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <assert.h>
#include <stdio.h>

// Aliases for types defined by the GNU Scientific Library
typedef gsl_vector vector_t;
typedef gsl_matrix matrix_t;


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
exit_t gsl_matrix_flatten (const gsl_matrix* mat, gsl_vector** vec);


/**
 * Print the contents of a gsl_vector to stdout.
 *
 * Parameters:
 *     vec - the vector to print
 *
 * Returns:
 *     void
 */
void gsl_vector_printf (const gsl_vector* vec);


/**
 * Print the conents of a gsl_matrix to stdout.
 *
 * Parameters:
 *     mat - the matrix to print
 *
 * Returns:
 *     void
 */
void gsl_matrix_printf (const gsl_matrix* mat);


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
exit_t gsl_matrix_get_rows (const gsl_matrix *mat, const gsl_vector *indices, gsl_matrix **result);


#endif
