#ifndef _CELLDATA_C_
#define _CELLDATA_C_

#include "celldata.h"


celldata* celldata_alloc () {
  return (celldata*) malloc (sizeof(celldata));
}

void celldata_init (celldata** data, const size_t ncells, const size_t ncoords, const size_t ndims) {
  
  size_t nrows = ncells * ncoords;
  size_t ncols = ndims;

  printf ("~> Creating celldata with <ncells=%lu ncoords=%lu ndims=%lu matrix %lu X %lu>\n", ncells, ncoords, ndims, nrows, ncols);
  *data = celldata_alloc();
  (*data)->ncells	= ncells;
  (*data)->ncoords = ncoords;
  (*data)->ndims   = ndims;
  (*data)->matrix  = gsl_matrix_calloc(ncells*ncoords, ndims);
}




matindex celldata_matindex (const celldata* cells, const cindex cell) {
  size_t i = cells->ncells+1 * cell;
  if (i >= cells->matrix->size1) {
    printf ("~> FAILURE: cell = %lu ncells = %lu i = %lu >= %lu\n", cell, cells->ncells, i,  cells->matrix->size1);
    celldata_printinfo ( cells );
  }
  assert (i < cells->matrix->size1);
  return i;
}

cindex celldata_cell_cindex (const celldata* cells, const matindex ix) {
  return ix % cells->ncells;
}


gsl_matrix celldata_get_cell (const celldata* data, const cindex cell) {
  size_t \
    row = celldata_matindex (data, cell),
    col = 0,
    nrs = data->ncoords,
    ncs = data->ndims;
  /* printf ("~> celldata_get_cell: slice  <row=%lu col=%lu nrows=%lu ncols=%lu>\n", row, col, nrs, ncs); */
  /* printf ("~> celldata_get_cell: actual <%lu X %lu>\n", data->matrix->size1, data->matrix->size2); */
  return gsl_matrix_submatrix (data->matrix, row, col, nrs, ncs).matrix;
}

gsl_vector celldata_get_coords (const celldata* data, const cindex cell, const cindex coord) {
  gsl_matrix mat = celldata_get_cell (data, cell);
  return gsl_matrix_row (&mat, coord).vector;
}

double celldata_get_value (const celldata* data, const cindex cell, const cindex coord, const cindex dim) {
  gsl_matrix mat = celldata_get_cell (data, cell);
  return gsl_matrix_get (&mat, coord, dim);
}

void celldata_set_value (celldata* data, const cindex cell, const cindex coord, const cindex dim, const double value) {
  matindex r = celldata_matindex (data, cell);
  r += coord;
  gsl_matrix_set (data->matrix, r, dim, value);
}
  
  
void celldata_printinfo(const celldata* cells) {
  printf ("celldata_printinfo\n");
  printf ("<celldata: ncells = %lu ncoords = %lu ndims = %lu matrix dimensions = %lu X %lu>\n",
	  cells->ncells, cells->ncoords, cells->ndims, cells->matrix->size1, cells->matrix->size2);
}

void celldata_printf (const celldata* cells) {
  celldata_printinfo (cells);
  for (int c=0; c<cells->ncells; c++){
    for (int x=0; x<cells->ncoords; x++){
      for (int d=0; d<cells->ndims; d++){
  	printf ("  [ %5d, %5d, %5d ] = %9.3f\n", c, x, d, celldata_get_value (cells, c, x, d)) ;
      }}}
}



void test_celldata_init () {
  celldata* cells;
  size_t ncells = 5;
  size_t ncoords = 9;
  size_t ndims = 3;

  printf ("test 1 %d\n", cells == NULL);

  celldata_init ( &cells, ncells, ncoords, ndims );
  printf ("test 2 %d\n", cells == NULL);
  celldata_printinfo ( cells );

  int counter = 0;


}



#endif
