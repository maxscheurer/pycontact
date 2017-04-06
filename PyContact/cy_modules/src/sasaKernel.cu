#include <stdio.h>

void __global__ SasaKernel(int natoms, float pairdist, const int npts,
    const int neighboursInSmem, float srad, float3* pos,
    float* radius, float3* points, float* sasa)
{
    extern __shared__ float4 s_mem[];
    const int thread = blockDim.x * blockIdx.x + threadIdx.x;
    for (int i = 0; i < neighboursInSmem; i++) {
        float4 atom;
        atom.x = pos[i].x;
        atom.y = pos[i].y;
        atom.z = pos[i].z;
        atom.w = radius[i];
        s_mem[i] = atom;
    }
    __syncthreads();

    if (thread < natoms) {
        float3 center = pos[thread];
        float currentRadius = radius[thread];

        float3 scaledPoint;
        int remainingPoints = npts;

        for (int i = 0; i < npts; i++) {
            scaledPoint.x = points[i].x * (currentRadius + srad) + center.x;
            scaledPoint.y = points[i].y * (currentRadius + srad) + center.y;
            scaledPoint.z = points[i].z * (currentRadius + srad) + center.z;

            bool overlap = false;

            for (int atomId = 0; atomId < natoms; atomId++) {
                if (atomId != thread) {
                    if (atomId < neighboursInSmem) {
                        float4 neighbour = s_mem[atomId];
                        float neighbourRad2 = neighbour.w + srad;
                        neighbourRad2 *= neighbourRad2;
                        float3 dr;
                        dr.x = scaledPoint.x - neighbour.x;
                        dr.y = scaledPoint.y - neighbour.y;
                        dr.z = scaledPoint.z - neighbour.z;

                        if (dr.x*dr.x + dr.y*dr.y + dr.z*dr.z <= neighbourRad2) {
                            overlap = true;
                            break;
                        }
                    } else {
                        float neighbourRad2 = radius[atomId] + srad;
                        neighbourRad2 *= neighbourRad2;
                        float3 neighbourCenter = pos[atomId];
                        float3 dr;
                        dr.x = scaledPoint.x - neighbourCenter.x;
                        dr.y = scaledPoint.y - neighbourCenter.y;
                        dr.z = scaledPoint.z - neighbourCenter.z;

                        if (dr.x*dr.x + dr.y*dr.y + dr.z*dr.z <= neighbourRad2) {
                            overlap = true;
                            break;
                        }
                    }
                }
            }
            if (overlap) {
                remainingPoints--;
            }
        }
        sasa[thread] = 12.5663706144 * powf(currentRadius+srad, 2) * (float)(remainingPoints) / npts;
    }
}
