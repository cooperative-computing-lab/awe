#ifndef _ASSIGN_H_
#define _ASSIGN_H_

#include "celldata.h"
#include "xdr_util.h"
#include "rmsd_calc.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>



int load_data (const char* path, celldata* data);




#endif
