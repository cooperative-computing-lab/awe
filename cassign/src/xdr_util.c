/**
 * xdr_util.c
 * 
 * Functions for reading and filtering atom information in XTC
 * files and xdrfile structures.
 */

#include "xdr_util.h"

/**
 * Print the summary of an XDRFrame's contents to stdout, including the number
 * of atoms in the frame and the dimensions of the coordinates.
 *
 * Parameters:
 *     frame - an allocated or initialized xdrframe pointer
 *
 * Returns:
 *     void
 */
void xdrframe_printsummary(const xdrframe* frame) {
  printf("<xdrframe: natoms = %lu matrix = %lu X %lu>",
          frame->natoms, frame->coords->size1, frame->coords->size2);
}


/**
 * Print the coordinate matrix of an XDRFrame to stdout.
 *
 * Parameters:
 *     frame - an allocated or initialized xdrframe pointer
 *
 * Returns:
 *     void
 */
void xdrframe_printf(const xdrframe* frame) {
  xdrframe_printsummary(frame); // Print the summary as a header
  printf("{\n");

  // Iterate over each row (atom) in the matrix
  for (int i=0; i<frame->coords->size1; i++){
    printf("    [ ");
    
    // Iterate over each atomic coordinate
    for (int j=0; j<frame->coords->size2; j++){
      printf("%.3f ", gsl_matrix_get (frame->coords, i, j));
    }

    printf("]\n");
  }

  printf("}\n");
}


/**
 * Allocate space for an xdrframe struct pointer.
 *
 * Parameters:
 *     None
 *
 * Returns:
 *     An XDRFrame pointer that has been allocated space on the heap
 */
xdrframe* xdrframe_alloc() {
  return (xdrframe*) malloc(sizeof(xdrframe));
}


/**
 * Initialize an XDRFrame pointer on the heap with a specified number of atoms.
 *
 * Parameters:
 *     frame  - a pointer to an xdrframe pointer where the initialized xdrframe
 *              pointer will be stored
 *     natoms - the number of rows to allocate to the coordinate matrix
 */
void xdrframe_init(xdrframe** frame, const size_t natoms) {
  *frame = xdrframe_alloc ();
  (*frame)->natoms = natoms;
  (*frame)->coords = gsl_matrix_calloc (natoms, XDR_DIM);
}


/**
 * Add coordinates to an initialized XDRFrame.
 *
 * Parameters:
 *     x     - an array of coordinates to set in the coordinate matrix
 *     frame - an initialized xdrframe pointer in which to set the coordinates
 */
void xdrframe_from_rvec(const rvec* x, xdrframe* frame) {
  // Ensure that the frame is of the correct size to hold the coordinates
  assert(frame->natoms == frame->coords->size1);
  assert(frame->coords->size2 == XDR_DIM);

  // Add each coordinate to the coortdinate matrix
  for (int i=0; i<frame->natoms; i++) {
    for (int j=0; j<XDR_DIM; j++) {
      gsl_matrix_set(frame->coords, i, j, x[i][j]);
    }
  }
}


/**
 * Get the last frame of an XTC trajectory file.
 *
 * Parameters:
 *     filename - filepath to a GROMACS XTC file to read
 *     frame    - the xdrframe in which to store the last from in the XTC file
 *
 * Returns:
 *     The answer to life, the universe, and everything or an error flag
 */
int xdrframe_last_in_xtc (const char* filename, xdrframe** frame) {

  printf("~> Loading last frame from xtc file: %s\n", filename);

  // Open the filepath and ensure it exists
  XDRFILE* file = xdrfile_open(filename, "r");
  if (file == NULL) {
    char emsg[50];
    sprintf(emsg, "~> Error opening xtc file: %s", filename);
    perror(emsg);
    return exdrFILENOTFOUND;
  }

  // Determine how many atoms are in the molecule and 
  // return if a read error occurs.
  int natoms;
  int result;
  if ( (result = read_xtc_natoms((char*)filename, &natoms)) != exdrOK ) {
    printf("~> Error reading number of atoms from: %s\n", filename);
    return result;
  }

  printf("~> Number of atoms: %d\n", natoms);

  // Setup variables needed to find the last frame.
  int step;
  float time;
  xdr_matrix box;
  float prec    = XTC_PRECISION;
  int frameix   = 0;
  xdr_vec* x    = (xdr_vec*) malloc (natoms*sizeof(xdr_vec));

  // Read and store coordinates until the end of the file is reached
  while ( (result = read_xtc(file, natoms, &step, &time, box, x, &prec)) == exdrOK ) {
    frameix++;
  }

  xdrfile_close (file);
  printf ("~> Read %d frames from: %s\n", frameix, filename);  

  // Allocate, initialize, and populate the frame.
  *frame = xdrframe_alloc ();
  xdrframe_init (frame, natoms);
  xdrframe_from_rvec (x, *frame);

  // Enlightenment has been reached. Let everyone know.
  return 42;
}


/**
 * Get a set of atom coordinates from a populated XDRFrame
 *
 * Parameters:
 *     frame    - a populated xdrframe pointer from which to read coordinates
 *     indices  - the list of atoms for which coordinates are needed
 *     newframe - a pointer to the xdrframe pointer in which to store the
 *                coordinates
 *
 * Returns:
 *     A flag indicating the read was successful
 */
exit_t xdrframe_select_atoms (const xdrframe *frame, const gsl_vector *indices, xdrframe **newframe) {
  // Ensure that the number of atoms to read is not larger than the number
  // of atoms in the molecule.
  assert (indices->size <= frame->natoms);

  // Initialize a new xdrframe and populate it with specified the coordinates
  xdrframe_init (newframe, indices->size);
  gsl_matrix_get_rows (frame->coords, indices, &(*newframe)->coords);
  return exitOK;
}
