
#include "exit_codes.h"


const char* exit_message[1 + exitSANITY] = {
  "OK",
  "Index out of bounds",
  "Bad header section for cells data file",
  "Bad data section for cells data file",
  "Sanity check failed"
};
