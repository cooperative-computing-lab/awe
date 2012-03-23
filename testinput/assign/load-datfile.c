
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>



void my_sscanf(){
  char string[] = "2.045 0.706 -0.403";
  float x, y, z;

  sscanf (string, "%f %f %f", &x, &y, &z);
  printf ("loaded x=%f y=%f z=%f\n", x, y, z);

}


int main (void) {

  my_sscanf();

  return 0;
}
