#ifndef _CELLDATA_H_
#define _CELLDATA_H_

#include "exit_codes.h"
#include "gsl_util.h"

#include <gsl/gsl_matrix.h>

#include <assert.h>
#include <stdlib.h>

typedef gsl_matrix cell_t;

typedef struct {
  size_t ncells, ncoords, ndims;
  cell_t *cells;
} celldata;

cell_t cell_t_malloc (const size_t ncoords, const size_t ndims);

void celldata_init (celldata** data, const size_t ncells, const size_t ncoords, const size_t ndims);

cell_t celldata_get_cell (const celldata* data, const size_t cell);

exit_t celdlata_set_cell (celldata *celldata, const size_t c, const cell_t cell);

vector_t celldata_get_coords (const celldata* data, const size_t cell, const size_t coord);

double celldata_get_value (const celldata* data, const size_t cell, const size_t coord, const size_t dim);

void celldata_set_value (celldata* data, const size_t cell, const size_t coord, const size_t dim, const double value);

void celldata_printinfo (const celldata* cells);

void celldata_printf (const celldata* cells);

exit_t celldata_load_file (const char* path, celldata** data);

exit_t celldata_get_rows (const celldata *cells, const vector_t *atomindices, celldata **newcells);

#endif
