/**
 * celldata.h
 *
 * Utilities for creating and managing celldata structures.
 */

#ifndef _CELLDATA_H_
#define _CELLDATA_H_

#include "exit_codes.h"
#include "gsl_util.h"

#include <gsl/gsl_matrix.h>

#include <assert.h>
#include <stdlib.h>

// An alias of gsl_matrix to distinguish code is dealing with cells
typedef gsl_matrix cell_t;

// A structure for dealing with all cells in a clustering model
typedef struct {
  size_t ncells, ncoords, ndims;
  cell_t *cells;
} celldata;


/**
 * Allocate space for a celldata struct on the heap.
 *
 * Parameters:
 *     None
 *
 * Returns:
 *     A pointer to a celldata structure
 */
celldata* celldata_alloc ();


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
void celldata_init (celldata** data, const size_t ncells, const size_t ncoords, const size_t ndims);


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
cell_t celldata_get_cell (const celldata* data, const size_t cell);


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
exit_t celdlata_set_cell (celldata *celldata, const size_t c, const cell_t cell);


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
vector_t celldata_get_coords (const celldata* data, const size_t cell, const size_t coord);


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
double celldata_get_value (const celldata* data, const size_t cell, const size_t coord, const size_t dim);


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
void celldata_set_value (celldata* data, const size_t cell, const size_t coord, const size_t dim, const double value);


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
void celldata_printinfo (const celldata* cells);


/**
 * Print out all information stored in a celldata structure to stdout.
 *
 * Parameters:
 *     celldata - the structure about which to print information
 *
 * Returns:
 *     void
 */
void celldata_printf (const celldata* cells);


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
exit_t celldata_load_file (const char* path, celldata** data);


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
exit_t celldata_get_rows (const celldata *cells, const vector_t *atomindices, celldata **newcells);


#endif
