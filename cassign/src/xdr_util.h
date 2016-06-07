/**
 * xdr_util.h
 * 
 * Functions and definitions for reading and filtering atom information in XTC
 * files and xdrfile structures.
 */

#ifndef _XDR_UTIL_H_
#define _XDR_UTIL_H_

#include "gsl_util.h"

#include "common.h"

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <gsl/gsl_matrix.h>

#include <stdlib.h>
#include <assert.h>

#define XTC_PRECISION 1000
#define XDR_DIM DIM

// Aliases for types defined by the GROMACS XDRFile library
typedef matrix xdr_matrix;
typedef rvec   xdr_vec;

// A structure for storing atom coordinate information
typedef struct {
  size_t natoms;      // The number of rows (atoms) in the matrix
  gsl_matrix* coords;
} xdrframe;


// See xdr_util.c for  more detailed descriptions

/**
 * Allocate space for an xdrframe struct pointer.
 */
xdrframe* xdrframe_alloc (void);

/**
 * Initialize an XDRFrame pointer on the heap with a specified number of atoms.
 */
void xdrframe_init (xdrframe** frame, const size_t natoms);

/**
 * Add coordinates to an initialized XDRFrame.
 */
void xdrframe_from_rvec (const rvec* x, xdrframe* frame);

/**
 * Get the last frame of an XTC trajectory file.
 */
int xdrframe_last_in_xtc (const char* filename, xdrframe** frame);

/**
 * Get a set of atom coordinates from a populated XDRFrame
 */
exit_t xdrframe_select_atoms (const xdrframe *frame, const gsl_vector *indices, xdrframe **newframe);

/**
 * Print the summary of an XDRFrame's contents to stdout, including the number
 * of atoms in the frame and the dimensions of the coordinates.
 */
void xdrframe_printsummary (const xdrframe* frame);

/**
 * Print the coordinate matrix of an XDRFrame to stdout.
 */
void xdrframe_printf (const xdrframe* frame);


#endif
