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

exit_t celldata_load_file (const char* path, celldata** data) {

  int BUFFER_SIZE = 100;
  char buffer [BUFFER_SIZE];
  int ncells, ncoords, ndims;

  FILE * file = fopen(path, "r");
  if (file == NULL) perror ("Error opening file");
  else {

    printf ("Reading header\n");

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of cells");
    sscanf (buffer, "ncells: %d", &ncells);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of coordinates");
    sscanf (buffer, "ncoords: %d", &ncoords);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of dimensions");
    sscanf (buffer, "ndims: %d", &ndims);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL) perror ("Error clearing empty line");

    printf("ncells = %d ncoords = %d ndims = %d\n", ncells, ncoords, ndims);
    celldata_init (data, ncells, ncoords, ndims);



    int lineno  = 0;
    int cell	= 0;
    int coord	= 0;
    int dim	= 0;
    float val	= 0;

    /* read in the data */
    printf ("Reading data...\n");
    while ( ! feof (file) ) {
      if ( fgets (buffer, BUFFER_SIZE, file) != NULL ){
	lineno += 1;
	if ( sscanf (buffer, "%f", &val) != 1 ) {
	  printf ("Failed to parse line %d of data section: '%s'\n", lineno, buffer);
	  exit (EXIT_FAILURE);
	}
	celldata_set_value (data, cell, coord, dim, val);

	/* printf ("line %d data[%d][%d][%d] = %f\n", lineno, cell, coord, dim, val); */

	/* printf ("incrementing dim %d\n", dim); */
	dim += 1;

	if (dim >= ndims){
	  /* printf ("reseting dim %d to 0\n", dim); */
	  dim = 0;
	}

	if (dim == 0) {
	  /* printf ("incrementing coord %d\n", coord); */
	  coord += 1;
	}

	if (coord >= ncoords) {
	  /* printf ("reseting coord %d\n", coord); */
	  coord = 0;
	}

	if (coord == 0 && dim == 0) {
	  /* printf ("incrementing cell %d\n", cell); */
	  cell += 1;
	}

	if (cell >= ncells) {
	  printf ("loaded %d cells\n", ncells);
	  break;
	}

      }
    }

    fclose(file);
  }
  printf ("...done\n");

  return exitOK;
}

#endif
