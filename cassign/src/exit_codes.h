/**
 * exit_codes.h
 * 
 * Utilities for determining function return codes for debugging and integrity.
 */

#ifndef _EXIT_CODES_H_
#define _EXIT_CODES_H_

// An enumeration of exit messages used as an index for exit_message (below)
typedef enum { 
  exitOK,
  exitPATH_NOT_FOUND,
  exitFAILURE,
  exitIX_OUT_OF_BOUNDS,
  exitCELLS_HEADER,
  exitCELLS_DATA,
  exitSANITY
} exit_t;

// An array of exit status messages indexed by exit_t
const char* exit_message[1 + exitSANITY];


#endif
