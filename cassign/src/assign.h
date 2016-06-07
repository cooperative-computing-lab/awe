/**
 * assign.h
 *
 * The main driver for the AWE-Assign program. Assigns the end of an MD
 * trajectory to a specific cluster in the conformation space of the
 * molecule.
 */

#ifndef _ASSIGN_H_
#define _ASSIGN_H_

#include "celldata.h"
#include "xdr_util.h"
#include "rmsd_calc.h"
#include "gsl_util.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>


/**
 * Load the indices of matrix rows to evaluate in the matrix of atomic
 * coordinates of a molecular structure.
 *
 * Parameters:
 *     mndxpath - the path to the file containing indices
 *     target   - the vector that contains the indices to examine
 *
 * Returns:
 *     An exit status representing the success or failure of the operation
 */
exit_t load_atomindices (const char *mndxpath, gsl_vector **target);


#endif
