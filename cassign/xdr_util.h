#ifndef _XDR_UTIL_H_
#define _XDR_UTIL_H_

#include "gsl_util.h"

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <gsl/gsl_matrix.h>

#include <stdlib.h>
#include <assert.h>

#define XTC_PRECISION 1000
#define XDR_DIM DIM

typedef matrix xdr_matrix;
typedef rvec   xdr_vec;

typedef struct {
  size_t natoms;
  gsl_matrix* coords;
} xdrframe;


xdrframe* xdrframe_alloc (void);
void xdrframe_init (xdrframe** frame, const size_t natoms);

void xdrframe_from_rvec (const rvec* x, xdrframe* frame);

int xdrframe_last_in_xtc (const char* filename, xdrframe** frame);

exit_t xdrframe_select_atoms (const xdrframe *frame, const gsl_vector *indices, xdrframe **newframe);

void xdrframe_printsummary (const xdrframe* frame);
void xdrframe_printf (const xdrframe* frame);


#endif
