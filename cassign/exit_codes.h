#ifndef _EXIT_CODES_H_
#define _EXIT_CODES_H_

typedef enum { 
  exitOK, exitPATH_NOT_FOUND, exitFAILURE, exitIX_OUT_OF_BOUNDS, exitCELLS_HEADER, exitCELLS_DATA, exitSANITY
} exit_t;

const char* exit_message[1 + exitSANITY];


#endif
