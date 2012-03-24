#ifndef _RMSD_CALC_H_
#define _RMSD_CALC_H_

#include "exit_codes.h"
#include "xdr_util.h"
#include "celldata.h"
#include "gsl_util.h"

#include "theobald_rmsd.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

typedef struct {
  size_t nreal, npad, ndim;
  double g;
  gsl_matrix* coords;
} theodata;


theodata* theodata_alloc ();

exit_t theodata_init (const size_t nreal, const size_t ndim, theodata** theo);

void theodata_printf (const theodata* theo);


/* Prepare a structure for rmsd computation using Theobald metric
     mat: IN: the 2-rank matrix of coordinates
     theo: OUT: the prepared data
*/
exit_t prepare_data (const gsl_matrix* mat, theodata** theo);



/* Center the structure *in place*
   - mat: IN, OUT: the 2-rank matrix of coordinates
*/
exit_t center_structure (gsl_matrix* mat);

exit_t flatten (const gsl_matrix* mat, gsl_vector* vec);

double compute_rmsd (const gsl_matrix* m1, const gsl_matrix* m2);

/* Compute the G value for input into TheoRMSD */
double calculate_theo_g (const gsl_matrix* mat);

#endif
