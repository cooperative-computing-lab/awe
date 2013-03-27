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



exit_t load_atomindices (const char *mndxpath, gsl_vector **target);




#endif
