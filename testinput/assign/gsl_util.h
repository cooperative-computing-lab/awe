#ifndef _GSL_UTIL_H_
#define _GSL_UTIL_H_

#include "exit_codes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <stdio.h>

exit_t gsl_matrix_flatten (const gsl_matrix* mat, gsl_vector** vec);

void gsl_vector_printf (const gsl_vector* vec);

void gsl_matrix_printf (const gsl_matrix* mat);

#endif
