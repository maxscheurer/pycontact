/*
This is the central piece of code. This file implements a class
(interface in gpuadder.hh) that takes data in on the cpu side, copies
it to the gpu, and exposes functions (increment and retreive) that let
you perform actions with the GPU

This class will get translated into python via swig
*/

#include <sasaKernel.cu>
#include <sasaCuda.hh>
#include <assert.h>
#include <iostream>


double calculate_sasa_cuda(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs,
                           const float *radius,const int npts, double srad, int pointstyle, int restricted,
                           const int* restrictedList) {
    std::cout << "Cuda Test!" << std::endl;
    NullKernel<<< 1, 1024 >>>();
    return 0.0f;
}
