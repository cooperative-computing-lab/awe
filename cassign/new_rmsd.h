#include "gsl_util.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <assert.h>
#include <math.h>


double rmsd (const gsl_matrix *A, const gsl_matrix *B);
