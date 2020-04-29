
#include <stdio.h>
#include <stdlib.h>
#include <cutil.h>

extern "C" void getlasterror_ (int *tag) {

  cudaError_t status;
  status = cudaGetLastError();
  printf("Last error: %d %d \n",*tag,status);
  printf("Cuda error: %s \n", cudaGetErrorString( status) );
  fflush(stdout);
  return;
}
