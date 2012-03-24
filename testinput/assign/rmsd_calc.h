#ifndef _RMSD_CALC_H_
#define _RMSD_CALC_H_

#include "exit_codes.h"
#include "xdr_util.h"
#include "celldata.h"

#include "theobald_rmsd.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>



/* Center the structure *in place*
   - mat: IN, OUT: the 2-rank matrix of coordinates
*/
exit_t center_structure (gsl_matrix* mat);

#endif
