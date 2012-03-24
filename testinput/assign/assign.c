
#include "assign.h"


int load_data (const char* path, celldata* data) {

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
    celldata_init (&data, ncells, ncoords, ndims);



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

 
  return 0;
}




/* float compute_rmsd (float** ref, float** structure, int natoms, int ndims) { */

/*   float *AData, *BData; */
/*   return 42; */

/*   /\* float msd = ls_rmsd2_aligned_T_g(nrealatoms,npaddedatoms,rowstride,(AData+i*truestride),BData,GAData[i],G_y); *\/ */

/* } */
  


int main (void) {

  const char* cells_file = "Gens.dat";
  const char* xtc_file   = "traj.xtc";

  celldata* cell_data;
  xdrframe* frame;

  load_data (cells_file, cell_data);
  xdrframe_last_in_xtc (xtc_file, &frame);
  xdrframe_printsummary (frame);
  xdrframe_printf (frame);

  center_structure (frame->coords);
  printf ("~> centered conformation\n");
  xdrframe_printsummary (frame);
  xdrframe_printf (frame);
    

  exit (EXIT_SUCCESS);
}
