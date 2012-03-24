#ifndef _CELLDATA_H_
#define _CELLDATA_H_

#include "exit_codes.h"

#include <gsl/gsl_matrix.h>

#include <assert.h>
#include <stdlib.h>

typedef struct {
  size_t ncells, ncoords, ndims;
  gsl_matrix* matrix;
} celldata;

typedef size_t cindex;
typedef size_t matindex;


celldata* celldata_alloc (void);

void celldata_init (celldata** data, const size_t ncells, const size_t ncoords, const size_t ndims);

matindex celldata_matindex (const celldata* cells, const cindex cell);

cindex celldata_cell_cindex (const celldata* cells, const matindex ix);

gsl_matrix celldata_get_cell (const celldata* data, const cindex cell);

gsl_vector celldata_get_coords (const celldata* data, const cindex cell, const cindex coord);

double celldata_get_value (const celldata* data, const cindex cell, const cindex coord, const cindex dim);

void celldata_set_value (celldata* data, const cindex cell, const cindex coord, const cindex dim, const double value);

void celldata_printinfo (const celldata* cells);

void celldata_printf (const celldata* cells);

exit_t celldata_load_file (const char* path, celldata* data);

#endif
