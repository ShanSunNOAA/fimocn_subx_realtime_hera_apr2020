// Routine to initialize the GPU
// Author:  Jacques Middlecoff
// Date:  September 2010 
// For Fortran this routine does nothing except return error=0.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cutil.h>
#include "ftocmacros.h"

extern "C" void gpuinit_ (int *npes,int *me,int *error) {

int argc=2;
char *argv0[]= {"","-device=0"};
char *argv1[]= {"","-device=1"};
cudaDeviceProp deviceProp;

*error = 0;

#if CUDART_VERSION < 2020
#error "This CUDART version does not support mapped memory!\n"
#endif

// Get properties and verify device 0 supports mapped memory
*error = cudaGetDeviceProperties(&deviceProp, 0);
if(*error != cudaSuccess) {
  printf("GPUinit.cu: cudaGetDeviceProperties error %d \n",*error);
  return;
}
if(!deviceProp.canMapHostMemory) {
  printf("GPUinit.cu: Device %d cannot map host memory!\n", 0);
  *error = -88;
  return;
}

if(*me%2 == 0)
{
  CUT_DEVICE_INIT(argc, argv0);
  printf("Processor %d %s \n",*me,argv0[1]);
} 
else
{
  CUT_DEVICE_INIT(argc, argv1);
  printf("Processor %d %s \n",*me,argv1[1]);
}

// set the device flags for mapping host memory
*error = cudaSetDeviceFlags(cudaDeviceMapHost);
if(*error != cudaSuccess) {
  printf("GPUinit.cu: cudaSetDeviceFlags error %d \n",*error);
  return;
}
}
