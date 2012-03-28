#ifndef _CELLDATA_C_
#define _CELLDATA_C_

#include "celldata.h"


celldata* celldata_alloc () {
  return (celldata*) malloc (sizeof(celldata));
}

void celldata_init (celldata** data, const size_t ncells, const size_t ncoords, const size_t ndims) {
  
  printf ("~> Creating celldata with <ncells=%lu ncoords=%lu ndims=%lu>\n", ncells, ncoords, ndims);
  *data = celldata_alloc();
  (*data)->ncells  = ncells;
  (*data)->ncoords = ncoords;
  (*data)->ndims   = ndims;
  (*data)->cells   = (cell_t*) malloc (ncells*sizeof(cell_t));

  for (int c=0; c<ncells; c++) {
    (*data)->cells[c] = * gsl_matrix_calloc (ncoords, ndims);
  }

}

cell_t celldata_get_cell (const celldata* data, const size_t cell) {
  assert (cell < data->ncells);
  return data->cells[cell];
}

vector_t celldata_get_coords (const celldata* data, const size_t cell, const size_t coord) {
  cell_t mat = celldata_get_cell (data, cell);
  return gsl_matrix_row (&mat, coord).vector;
}

double celldata_get_value (const celldata* data, const size_t cell, const size_t coord, const size_t dim) {
  cell_t *mat = &data->cells[cell];
  return gsl_matrix_get (mat, coord, dim);
}

void celldata_set_value (celldata* data, const size_t cell, const size_t coord, const size_t dim, const double value) {
  cell_t *c = &data->cells[cell];
  gsl_matrix_set (c, coord, dim, value);
}

void celldata_printinfo(const celldata* cells) {
  printf ("<celldata: ncells = %lu ncoords = %lu ndims = %lu>", cells->ncells, cells->ncoords, cells->ndims);
}

void celldata_printf (const celldata* cells) {
  celldata_printinfo (cells);
  for (int i=0; i<cells->ncells; i++) {
    cell_t cell = celldata_get_cell (cells, i);
    gsl_matrix_printf (&cell);
  }
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
	celldata_set_value (*data, cell, coord, dim, (double) val);

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

exit_t celldata_set_cell (celldata *celldata, const size_t c, const cell_t *cell) {
  assert (c < celldata->ncells);
  assert (cell->size1 == celldata->ncoords);
  assert (cell->size2 == celldata->ndims);

  for (int r=0; r<celldata->ncoords; r++) {
    for (int col=0; col<celldata->ndims; col++) {
      const double v = gsl_matrix_get (cell, r, col);
      celldata_set_value (celldata, c, r, col, v);
    }}

  return exitOK;
}

exit_t celldata_get_rows (const celldata *cells, const vector_t *atomindices, celldata **newcells) {
  assert (atomindices->size < cells->ncoords);
  celldata_init (newcells, cells->ncells, atomindices->size, cells->ndims);

  for (int cid=0; cid<cells->ncells; cid++) {
    const cell_t oldcell = celldata_get_cell (cells, cid);
    cell_t *newcell;
    gsl_matrix_get_rows (&oldcell, atomindices, &newcell);
    celldata_set_cell (*newcells, cid, newcell);
  }

  return exitOK;
}  

#endif
