/**
 * exit_codes.c
 * 
 * Utilities for determining function return codes for debugging and integrity.
 */

#include "exit_codes.h"

/**
 * A mapping for the exit_t enumerated type to dynamically return function
 * results feedback. See exit_codes.h for the exit_t definition.
 */
const char* exit_message[1 + exitSANITY] = {
  "OK",
  "Path not found",
  "Failure",
  "Index out of bounds",
  "Bad header section for cells data file",
  "Bad data section for cells data file",
  "Sanity check failed"
};
