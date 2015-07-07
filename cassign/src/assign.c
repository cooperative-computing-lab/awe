/**
 * assign.c
 *
 * The main driver for the AWE-Assign program. Assigns the end of an MD
 * trajectory to a specific cluster in the conformation space of the
 * molecule.
 */

#include "assign.h"

/**
 * Load the indices of matrix rows to evaluate in the matrix of atomic
 * coordinates of a molecular structure.
 *
 * Parameters:
 *     mndxpath - the path to the file containing indices
 *     target   - the vector that contains the indices to examine
 *
 * Returns:
 *     An exit status representing the success or failure of the operation
 */
exit_t load_atomindices (const char *mndxpath, gsl_vector **target) {
  /** 1) get the number of lines the file
      2) allocate the vector
      3) reread the file, inserting values into the vector

      This access the file twice, but unfortunately I'm not (yet)
      aware of any standard list datastructure for C, and am willing
      to accept the redundancy at the moment
  */

  // There are probably OS utilities in place to determine the number of
  // lines in a file--these can be conditionally included depending on the
  // compiling OS/architecture.

  exit_t status = exitOK;

  int linecount = 0;
  size_t ndx;

  // Ensure that the file exists and can be opened
  FILE *mndx = fopen (mndxpath, "r");
  
  if (mndx == NULL) {
    return exitPATH_NOT_FOUND;
  }

  // Read each newline to the end of the file to get the number of lines
  while ( !(feof(mndx)) ) {
    int count = fscanf(mndx, "%d", &ndx);
    
    if (count == 1) {
      linecount ++;
    }
  }

  fclose(mndx);

  if ( !(status == exitOK) ) {
    return status;
  }

  // Reopen the file
  // The file pointer could also be rewound to the beginning of the file
  mndx = fopen(mndxpath, "r");

  // Allocate the output vector assuming one index per line
  *target = gsl_vector_calloc(linecount);
  size_t i = 0;

  // Read each index and add it to the output vector
  while ( !(feof(mndx)) ) {
    int count = fscanf(mndx, "%d", &ndx);
    
    // Only add to the vector if a value was read
    if (count == 1) {
      gsl_vector_set(*target, i, (int) ndx);
      i++;
    }

    // That isn't how fscanf works
    else if ( count > 1 ) {
      perror("Atom indices file incorrectly formatted");
      status = exitFAILURE;
      break;
    }

    // If no values were read, assume EOF and end the loop
    else {
      break;
    }
  }

  fclose(mndx);

  return status;
}


/**
 * The main driver for the assignment algorithm.
 */
int main (int argc, char *argv[]) {

  // This seems unneccessary. Rewrite not to use the extra pointers.
  const char
    *cells_file    = argv[1], //"Gens.dat",
    *cell_ndx_file = argv[2], //"AtomIndices.dat",
    *xtc_file      = argv[3], //"traj.xtc",
    *xtc_ndx_file  = argv[4], //"AtomIndices.dat",
    *out_file      = argv[5]; //"cell2.dat";

  printf("~> Cells file: %s\n", cells_file);
  printf("~> Xtc file: %s\n", xtc_file);

  // Create the variables to be used in 
  celldata* cell_data;
  xdrframe* frame;
  gsl_vector *cell_indices, *xtc_indices;

  // Load the structure containing the cell information.
  // This should include every coordinate of every atom in every cell.
  celldata_load_file(cells_file, &cell_data);

  // Load the indices to evaluate
  load_atomindices(cell_ndx_file, &cell_indices);
  printf("~> Cell Atom Indices: ");
  gsl_vector_printf(cell_indices);
  printf("\n");

  printf("~> Loaded cells: ");
  celldata_printinfo(cell_data);
  printf("\n");

  // Load the last frame of the trajectory
  xdrframe_last_in_xtc(xtc_file, &frame);

  xdrframe_printsummary(frame);
  printf("\n");

  // REDUNDANT THIS SHOULD BE THE SAME AS THE OTHER INDICES FILE
  load_atomindices(xtc_ndx_file, &xtc_indices);
  printf("~> XTC Atom Indices: ");
  gsl_vector_printf(xtc_indices);
  printf("\n");

  // Cut down on the number of atoms for RMSD
  celldata *newcells;
  xdrframe *newframe;
  celldata_get_rows(cell_data, cell_indices, &newcells);
  xdrframe_select_atoms(frame, xtc_indices , &newframe);

  if (newcells->ncoords != newframe->natoms) {
    printf("Number of atoms mismatch after selecting atom indices: cell = %lu structure = %lu\n", newcells->ncoords, newframe->natoms);
    exit(1);
  }

  // THE ORIGINAL MATRICES ARE STILL ALIVE. THEY CAN AND SHOULD BE FREED

  printf("~> Using: ");
  celldata_printinfo(newcells);
  printf("\n");
  printf("~> Using: ");
  xdrframe_printsummary(newframe);
  printf("\n");
  

  printf("~> Computing rmsds...\n");
  
  // Set up variables for RMSD calculations
  int assignment = 0;
  double
    minrmsd = DBL_MAX,
    maxrmsd = -1;
  int asns = 0;
  
  for (int c=0; c<cell_data->ncells; c++) {
    // INITIALIZING NEW VARS EACH TIME IS WASTEFUL (POSSIBLY)
    const gsl_matrix cell = celldata_get_cell(newcells, c);
    const gsl_matrix *coords = newframe->coords;
    const double rmsd = compute_rmsd(&cell, coords);
    printf("~~> rmsd to cell %3d: %8.5f\n", c, rmsd);
    
    // ARE THE MIN AND MAX NECESSARY?
    // GET RID OF MAX--IT'S UNUSED
    if (rmsd < minrmsd) {
      minrmsd = rmsd;
      assignment = c;
      asns++;
    }

    if (rmsd > maxrmsd ) {
      maxrmsd = rmsd;
    }
  }

  printf("~> Assignment: %d jumps: %d\n", assignment, asns);
  printf("~> minrmsd: %f maxrmsd: %f\n", minrmsd, maxrmsd);

  printf("~> Saving assignment to: %s\n", out_file);
  
  // Write the assignment to file if write is possible
  FILE *out = fopen(out_file, "w");
  
  if (out == NULL) {
    char emsg[50];
    sprintf(emsg, "Could not open for writing: %s\n", out_file);
    perror(emsg);
    exit(EXIT_FAILURE);
  }

  fprintf(out, "%d", assignment);
  fclose(out);

  exit(EXIT_SUCCESS);
}
