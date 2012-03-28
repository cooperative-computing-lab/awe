

#include "xdr_util.h"

void xdrframe_printsummary (const xdrframe* frame) {
  printf ("~> <xdrframe: natoms = %lu matrix = %lu X %lu>", frame->natoms, frame->coords->size1, frame->coords->size2);
}

void xdrframe_printf (const xdrframe* frame) {
  xdrframe_printsummary (frame);
  printf ("{\n");
  for (int i=0; i<frame->coords->size1; i++){
    printf ("    [ ");
    for (int j=0; j<frame->coords->size2; j++){
      printf ("%.3f ", gsl_matrix_get (frame->coords, i, j));
    }
    printf ("]\n");
  }
  printf ("}\n");
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


exit_t xdrframe_select_atoms (const xdrframe *frame, const gsl_vector *indices, xdrframe **newframe) {
  assert (indices->size <frame->natoms);

  xdrframe_init (newframe, indices->size);
  gsl_matrix_get_rows (frame->coords, indices, &(*newframe)->coords);
  return exitOK;
}


exit_t xdrframe_load_atomindices (const char *mndxpath, gsl_vector **target) {
  /** 1) get the number of lines the file
      2) allocate the vector
      3) reread the file

      This access the file twice, buf unfortunately, I'm not (yet)
      aware of any standard list datastructure for C, and am willing
      to accept the redundancy at the moment
  */

  exit_t status = exitOK;

  const int BUFFER_SIZE = 100;
  char buffer[BUFFER_SIZE];
  int linecount = 0;
  size_t ndx;

  FILE *mndx = fopen (mndxpath, "r");
  if (mndx == NULL) { return exitPATH_NOT_FOUND; }

  while ( ! feof (mndx) ) {
    int count = fscanf (mndx, "%d", &ndx);
    if (count == 1)
      { linecount ++; }
  }

  fclose (mndx);
  if ( ! status == exitOK ) { return status; }

  mndx = fopen (mndxpath, "r");

  *target = gsl_vector_calloc (linecount);
  size_t i = 0;

  while ( ! feof (mndx) ) {
    int count = fscanf (mndx, "%d", &ndx);
    if (count == 1) {
      gsl_vector_set (*target, i, (int)ndx);
      printf ("[i] = %f\n", gsl_vector_get (*target, i));
      i++;
    }
    else if ( count > 1 ) {
      perror ("Atom indices file incorrectly formatted");
      status = exitFAILURE;
      break;
    }
    else { break; }
  }

  fclose (mndx);

  return status;
}
