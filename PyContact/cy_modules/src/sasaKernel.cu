#include <stdio.h>

void __global__ SasaKernel(int natoms, float pairdist, const int npts, float srad, float3* pos,
    float* radius, float3* points, float* sasa)
{
    const int thread = blockDim.x * blockIdx.x + threadIdx.x;
    if (thread < natoms) {
        float3 center = pos[thread];
        float currentRadius = radius[thread];

        // Testing overlap for each point
        float3 scaledPoint;
        int remainingPoints = npts;

        for (int i = 0; i < npts; i++) {
            scaledPoint.x = points[i].x * (currentRadius + srad) + center.x;
            scaledPoint.y = points[i].y * (currentRadius + srad) + center.y;
            scaledPoint.z = points[i].z * (currentRadius + srad) + center.z;

            bool overlap = false;

            for (int atomId = 0; atomId < natoms; atomId++) {
                if (atomId != thread) {
                    float neighbourRad2 = radius[atomId] + srad;
                    neighbourRad2 *= neighbourRad2;
                    float3 neighbourCenter = pos[atomId];
                    float3 dr;
                    dr.x = scaledPoint.x - neighbourCenter.x;
                    dr.y = scaledPoint.y - neighbourCenter.y;
                    dr.z = scaledPoint.z - neighbourCenter.z;

                    // if (thread == 99 && atomId == 1000) {
                    //     printf("%f %f %f\n", dr.x, dr.y, dr.z);
                    // }

                    if (dr.x*dr.x + dr.y*dr.y + dr.z*dr.z <= neighbourRad2) {
                        overlap = true;
                        break;
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
