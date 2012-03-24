

#include "xdr_util.h"

void xdrframe_printsummary (const xdrframe* frame) {
  printf ("~> <xdrframe: natoms = %lu matrix = %lu X %lu>\n", frame->natoms, frame->coords->size1, frame->coords->size2);
}

void xdrframe_printf (const xdrframe* frame) {
  for (int i=0; i<frame->coords->size1; i++){
    printf ("  [ ");
    for (int j=0; j<frame->coords->size2; j++){
      printf ("%.3f ", gsl_matrix_get (frame->coords, i, j));
    }
    printf ("]\n");
  }
}


xdrframe* xdrframe_alloc () {
  return (xdrframe*) malloc (sizeof(xdrframe));
}

void xdrframe_init (xdrframe** frame, const size_t natoms) {
  *frame = xdrframe_alloc ();
  (*frame)->natoms = natoms;
  (*frame)->coords = gsl_matrix_calloc (natoms, XDR_DIM);
}


void xdrframe_from_rvec (const rvec* x, xdrframe* frame) {
  assert (frame->natoms == frame->coords->size1);
  assert (frame->coords->size2 == XDR_DIM);

  for (int i=0; i<frame->natoms; i++) {
    for (int j=0; j<XDR_DIM; j++) {
      gsl_matrix_set (frame->coords, i, j, x[i][j]);
    }}

}

int xdrframe_last_in_xtc (const char* filename, xdrframe** frame) {

  printf ("~> Loading last frame from xtc file: %s\n", filename);

  XDRFILE* file = xdrfile_open (filename, "r");
  if (file == NULL) {
    char emsg[50];
    sprintf (emsg, "~> Error opening xtc file: %s", filename);
    perror(emsg);
    return exdrFILENOTFOUND;
  }

  int natoms;
  int result;
  if ( (result = read_xtc_natoms ((char*)filename, &natoms) ) != exdrOK ) {
    printf ("~> Error reading number of atoms from: %s\n", filename);
    return result;
  }

  printf ("~> Number of atoms: %d\n", natoms);

  int step;
  float time;
  xdr_matrix box;
  float prec	= XTC_PRECISION;
  int frameix	= 0;
  xdr_vec* x	= (xdr_vec*) malloc (natoms*sizeof(xdr_vec));

  while ( (result = read_xtc (file, natoms, &step, &time, box, x, &prec)) == exdrOK )
    { frameix++; }
  xdrfile_close (file);
  printf ("~> Read %d frames from: %s\n", frameix, filename);

  *frame = xdrframe_alloc ();
  xdrframe_init (frame, natoms);
  xdrframe_from_rvec (x, *frame);

  return 42;

}

