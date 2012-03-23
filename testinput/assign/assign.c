
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>

#define XTC_PRECISION 1000
#define XDR_DIM       DIM


int load_data (const char* path, int* pncells, int* pncoords, int* pndims, float*** data) {

  int BUFFER_SIZE = 100;
  char buffer [BUFFER_SIZE];
  int ncells, ncoords, ndims;

  FILE * file = fopen(path, "r");
  if (file == NULL) perror ("Error opening file");
  else {

    printf ("Reading header\n");

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of cells");
    sscanf (buffer, "ncells: %d", pncells);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of coordinates");
    sscanf (buffer, "ncoords: %d", pncoords);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of dimensions");
    sscanf (buffer, "ndims: %d", pndims);
    ncells	= *pncells;
    ncoords	= *pncoords;
    ndims	= *pndims;

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL) perror ("Error clearing empty line");

    printf("ncells = %d ncoords = %d ndims = %d\n", ncells, ncoords, ndims);

    /* allocate and initialize the array */
    printf ("Allocating and initializing %d X %d X %d array\n", ncells, ncoords, ndims);
    data = (float***) malloc (ncells*sizeof(float**));
    for (int c=0; c<ncells; c++){
      data[c] = (float**) malloc (ncoords*sizeof(float*));
      for (int x=0; x<ncoords; x++){
	data[c][x] = (float*) malloc (ndims*sizeof(float));
	for (int d=0; d<ndims; d++){
	  data[c][x][d] = 0;
	}}}


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
	data[cell][coord][dim] = val;

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

  // sanity check
  for (int c=0; c<ncells; c++){
    for (int x=0; x<ncoords; x++){
      for (int d=0; d<ndims; d++){
	assert (data[c][x][d] != 0.0);
	assert (data[c][x][d] > 0 || data[c][x][d] < 0);
      }}}

  return 0;

}

int load_last_xtc_frame (const char* path, float*** coords, int* natoms) {

  printf ("Loading last frame from xtc file %s\n", path);

  XDRFILE * file = xdrfile_open(path, "r");
  if (file == NULL) {
    char emsg[100];
    sprintf (emsg, "Error opening xtc file %s", path);
    perror (emsg);
    return exdrFILENOTFOUND;
  }

  int result;

  if ( (result = read_xtc_natoms ((char*) path, natoms)) != exdrOK ) {
    printf ("Error reading number of atoms from %s", path);
    return result;
  }

  printf ("Number of atoms: %d\n", *natoms);

  int step;
  float time;
  matrix box;
  float prec	= 1000;
  int frame	= 0;
  rvec* x	= (rvec*) malloc ((*natoms)*sizeof(rvec));

  while ( (result = read_xtc (file, *natoms, &step, &time, box, x, &prec))  == exdrOK )
    { frame ++; }
  printf ("Read %d frames from %s\n", frame, path);

  (*coords) = (float**) malloc ((*natoms)*sizeof(float*));
  for (int a=0; a<(*natoms); a++){
    (*coords)[a] = (float*) malloc (XDR_DIM*sizeof(float));
    for (int d=0; d<XDR_DIM; d++) {
      (*coords)[a][d] = x[a][d];
    }
  }

  xdrfile_close(file);
  return exdrOK;


}

int main (void) {

  const char* cells_file = "Gens.dat";
  const char* xtc_file   = "traj.xtc";

  int ncells, ncoords, ndims;
  float*** data;
  load_data(cells_file, &ncells, &ncoords, &ndims, data);
  assert (ndims == 3);   // sanity check

  int natoms;
  float** coords;
  load_last_xtc_frame(xtc_file, &coords, &natoms);
  // sanity check  
  for (int a=0; a<natoms; a++){
    for (int d=0; d<XDR_DIM; d++) {
      assert (coords[a][d] != 0.0);
    }}

  exit (EXIT_SUCCESS);
}
