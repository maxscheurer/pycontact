#include <stdio.h>

void __global__ SasaKernel(int natoms, float pairdist, const int npts, float srad, float3* pos,
    float* radius, float3* points, float* sasa)
{
    sasa[5] = 42.0f;
}
