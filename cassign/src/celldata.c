/**
 * celldata.c
 *
 * Utilities for creating and managing celldata structures.
 */

#ifndef _CELLDATA_C_
#define _CELLDATA_C_

#include "celldata.h"


/**
 * Allocate space for a celldata struct on the heap.
 *
 * Parameters:
 *     None
 *
 * Returns:
 *     A pointer to a celldata structure
 */
celldata* celldata_alloc () {
  return (celldata*) malloc(sizeof(celldata));
}


/**
 * Initialize the contents of a celldata structure.
 *
 * Parameters:
 *     data    - a pointer to a celldata pointer in which the initialized
 *               structure will be stored
 *     ncells  - the number of cells to intialize 
 *     ncoords - the number of coordinates in each cell
 *     ndims   - the number of dimensions in each set of coordinates
 *
 * Returns:
 *     void
 */
void celldata_init (celldata** data, const size_t ncells, const size_t ncoords,
                    const size_t ndims) {
  
  printf("~> Creating celldata with <ncells=%lu ncoords=%lu ndims=%lu>\n",
         ncells, ncoords, ndims);
  
  // Allocate the celldata structure pointer
  *data = celldata_alloc();

  // Set the fields of the structure 
  (*data)->ncells  = ncells;
  (*data)->ncoords = ncoords;
  (*data)->ndims   = ndims;

  // Allocate space for the correct number of cells
  (*data)->cells   = (cell_t*) malloc(ncells*sizeof(cell_t));

  // Allocate an ncoords X ndims matrix for each cell
  for (int c=0; c<ncells; c++) {
    (*data)->cells[c] = * gsl_matrix_calloc(ncoords, ndims);
  }

}


/**
 * Get a cell from the celldata structure.
 *
 * Parameters:
 *     data - an initialized celldata structure from which to get a cell
 *     cell - the index of the cell to get from the celldata structure
 *
 * Returns:
 *     The cell_t instance at index cell
 */
cell_t celldata_get_cell (const celldata* data, const size_t cell) {
  // Ensure that the index a legal index
  assert(cell < data->ncells);
  return data->cells[cell];
}


/**
 * Get a specific set of coordinates from a cell.
 *
 * Parameters:
 *     data  - the celldata structure from which to get the coordinates
 *     cell  - the index of the cell from which to get coordinates
 *     coord - the index of the coordinates in the cell
 *
 * Returns:
 *     A vector_t containing the requested coordinates
 */
vector_t celldata_get_coords (const celldata* data, const size_t cell,
                              const size_t coord) {
  // Get the correct matrix (cell), then get the correct row (coordinates)
  cell_t mat = celldata_get_cell(data, cell);
  return gsl_matrix_row(&mat, coord).vector;
}


/**
 * Get one dimension from a particular set of coordinates from a particular
 * cell.
 *
 * Parameters:
 *     data  - the celldata structure from which to get the coordinates
 *     cell  - the index of the cell from which to get coordinates
 *     coord - the index of the coordinates in the cell
 *     dim   - the index of the dimension to get from the coordinate set
 *
 * Returns:
 *     The double at the specified cell coordinate set dimension.
 */
double celldata_get_value (const celldata* data, const size_t cell,
                           const size_t coord, const size_t dim) {
  // Get the correct matrix (cell), then get the dimension from the coordinates
  cell_t *mat = &data->cells[cell];
  return gsl_matrix_get(mat, coord, dim);
}


/**
 * Set one dimension of a particular set of coordinates of a particular cell.
 *
 * Parameters:
 *     data  - the celldata structure from which to get the coordinates
 *     cell  - the index of the cell from which to get coordinates
 *     coord - the index of the coordinates in the cell
 *     dim   - the index of the dimension to get from the coordinate set
 *     value - the new value of the dimension
 *
 * Returns:
 *     void
 */
void celldata_set_value (celldata* data, const size_t cell,
                         const size_t coord, const size_t dim,
                         const double value) {
  // Get the correct matrix (cell) and them set the intended dimension
  cell_t *c = &data->cells[cell];
  gsl_matrix_set(c, coord, dim, value);
}


/**
 * Print non-coordinate information about a celldata structure to stdout.
 * Prints the number of cells in the structure, number of coordinates in each
 * cell, and the dimensionality of each coordinate set.
 *
 * Parameters:
 *     celldata - the structure about which to print information
 *
 * Returns:
 *     void
 */
void celldata_printinfo(const celldata* cells) {
  printf("<celldata: ncells = %lu ncoords = %lu ndims = %lu>",
         cells->ncells, cells->ncoords, cells->ndims);
}


/**
 * Print out all information stored in a celldata structure to stdout.
 *
 * Parameters:
 *     celldata - the structure about which to print information
 *
 * Returns:
 *     void
 */
void celldata_printf (const celldata* cells) {
  // Print out the header info
  celldata_printinfo (cells);

  // Print out all information in all cells
  for (int i=0; i<cells->ncells; i++) {
    cell_t cell = celldata_get_cell(cells, i);
    gsl_matrix_printf(&cell);
  }
}


/**
 * Populate a celldata structure from the contents of a structured file.
 *
 * Parameters:
 *     path - the path to the file containing cell data
 *     data - the structure to populate
 * 
 * Returns:
 *     An exit status indicating the success of the operation.
 */
exit_t celldata_load_file (const char* path, celldata** data) {

  int BUFFER_SIZE = 100;
  char buffer [BUFFER_SIZE];
  int ncells, ncoords, ndims;

  // Open a the file containing cell definitions
  FILE * file = fopen(path, "r");
  
  // Ensure that the file was opened properly before proceeding
  if (file == NULL) {
    perror("Error opening file");
  }

  else {

    // The header should contain formatted information defining the shape of
    // the celldata->cells pointer.
    printf ("Reading header\n");

    // Get the number of cells
    if ( fgets(buffer, BUFFER_SIZE, file) == NULL ) {
      perror("Error reading number of cells");
    }

    sscanf(buffer, "ncells: %d", &ncells);

    // Get the number of coordinates per cell
    if ( fgets(buffer, BUFFER_SIZE, file) == NULL ) {
      perror("Error reading number of coordinates");
    }

    sscanf(buffer, "ncoords: %d", &ncoords);

    // Get the dimensionality of each coordinate
    if ( fgets(buffer, BUFFER_SIZE, file) == NULL ) {
      perror("Error reading number of dimensions");
    }

    sscanf(buffer, "ndims: %d", &ndims);

    // Move past the empty line after the shape information
    if ( fgets(buffer, BUFFER_SIZE, file) == NULL) {
      perror("Error clearing empty line");
    }

    // Prepare the structure for population
    printf("ncells = %d ncoords = %d ndims = %d\n", ncells, ncoords, ndims);
    celldata_init(data, ncells, ncoords, ndims);

    // Set up variables for reading in all of the dimension values
    int lineno  = 0; // Possibly a legacy variable. It doesn't seem to be used.
    int cell    = 0;
    int coord   = 0;
    int dim     = 0;
    float val   = 0;

    /* read in the data */
    // This is bad. This could (and probably should) be rewritten with loops.
    printf ("Reading data...\n");

    // Read until the end of the file
    while ( !(feof(file)) ) {

      // Get a line from the file
      if ( fgets(buffer, BUFFER_SIZE, file) != NULL ) {
        lineno += 1;

        // Get the contents of the line (a floating point value), or exit the
        // program if no value was read.
        if ( sscanf(buffer, "%f", &val) != 1 ) {
          printf("Failed to parse line %d of data section: '%s'\n",
                 lineno, buffer);
          exit(EXIT_FAILURE);
        }

        // Set the value at the proper index as a double
        celldata_set_value(*data, cell, coord, dim, (double) val);

        /* printf("line %d data[%d][%d][%d] = %f\n",
                  lineno, cell, coord, dim, val); */

        /* printf ("incrementing dim %d\n", dim); */

        // Increment dim and then determine how to update indices
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

        // Break out of the loop if there is more information but the
        // structure is full
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


/**
 * Set the state of a cell in a celldata structure.
 *
 * Parameter:
 *     celldata - the celldata structure containing the cell to update
 *     c        - the index of the cell to update
 *     cell     - the cell containing the new state information
 * 
 * Returns:
 *     An exit code representing the success of the flattening operation
 */
exit_t celldata_set_cell(celldata *celldata, const size_t c,
                         const cell_t *cell) {
  // Ensure that the cell exists and it is the same size as the update cell
  assert(c < celldata->ncells);
  assert(cell->size1 == celldata->ncoords);
  assert(cell->size2 == celldata->ndims);

  // Iteratively update the cell by row
  for (int r=0; r<celldata->ncoords; r++) {
    for (int col=0; col<celldata->ndims; col++) {
      const double v = gsl_matrix_get(cell, r, col);
      celldata_set_value(celldata, c, r, col, v);
    }
  }

  return exitOK;
}

/**
 * Filter cell definitions by certain indices and create a new celldata struct
 * containing the results of filtering.
 *
 * Parameter:
 *     cells       - the celldata structure containing the cell to get
 *                   information from
 *     atomindices - the indices of coordinates to get from each cell
 *     newcells    - the structure in which the resulting (filtered) cells will
 *                   be stored
 * 
 * Returns:
 *     An exit code representing the success of the operation
 */
exit_t celldata_get_rows(const celldata *cells, const vector_t *atomindices,
                         celldata **newcells) {
  // Ensure that the indices are in range
  assert(atomindices->size <= cells->ncoords);

  // Initialize the struct that will hold the result
  celldata_init(newcells, cells->ncells, atomindices->size, cells->ndims);

  // For each cell, filter by the rows in atomindices and add the filtered
  // cell to the resulting structure
  for (int cid=0; cid<cells->ncells; cid++) {
    const cell_t oldcell = celldata_get_cell(cells, cid);
    cell_t *newcell;
    gsl_matrix_get_rows(&oldcell, atomindices, &newcell);
    celldata_set_cell(*newcells, cid, newcell);
  }

  return exitOK;
}  

#endif
