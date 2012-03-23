
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#define BUFFER_SIZE 100


int main (void) {

  char buffer [BUFFER_SIZE];
  int ncells, ncoords, ndims;
  float ***data;

  FILE * file = fopen("Gens.dat", "r");
  if (file == NULL) perror ("Error opening file");
  else {

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of cells");
    sscanf (buffer, "ncells: %d", &ncells);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of coordinates");
    sscanf (buffer, "ncoords: %d", &ncoords);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL ) perror ("Error reading number of dimensions");
    sscanf (buffer, "ndims: %d", &ndims);

    if ( fgets (buffer, BUFFER_SIZE, file) == NULL) perror ("Error clearing empty line");

    printf("ncells = %d ncoords = %d ndims = %d\n", ncells, ncoords, ndims);

    /* allocate and initialize the array */
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

  for (int c=0; c<ncells; c++){
    printf("{%d", c);
    for (int x=0; x<ncoords; x++){
      printf("\t[ ");
      for (int d=0; d<ndims; d++){
	printf("%f ", data[c][x][d]);
      }
      printf("]\n");
    }
    printf("}\n");
  }


  exit (EXIT_SUCCESS);
}
