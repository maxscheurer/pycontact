#include <sasaKernel.cu>
#include <sasaCuda.hh>
#include <assert.h>
#include <iostream>


double calculate_sasa_cuda( float *pos, int natoms, float pairdist,  float *radius,  int npts,
                           float srad, int pointstyle, int restricted,  int* restrictedList) {

    const int blockSize = 256;
    const int gridSize = ceil(static_cast<float>(natoms) / static_cast<float>(blockSize));

    //std::cout << "Grid Size: " << gridSize << std::endl;
    //std::cout << "Block Size: " << blockSize << std::endl;

    float* h_sasa = NULL;
    float* d_sasa = NULL;
    float3* d_pos = NULL;
    float* d_radius = NULL;
    float3* h_points = NULL;
    float3* d_points = NULL;


    h_sasa = static_cast<float*>(malloc(static_cast<size_t>(natoms * sizeof(*h_sasa))));
    h_points = static_cast<float3*>(malloc(static_cast<size_t>(npts * sizeof(*h_points))));
    cudaMalloc(&d_sasa, static_cast<size_t>(natoms * sizeof(*d_sasa)));
    cudaMalloc(&d_pos, static_cast<size_t>(natoms * sizeof(*d_pos)));
    cudaMalloc(&d_radius, static_cast<size_t>(natoms * sizeof(*d_radius)));
    cudaMalloc(&d_points, static_cast<size_t>(npts * sizeof(*d_points)));

    memset(h_sasa, 0.0f, natoms * sizeof(*h_sasa));

    float phi_k = 0.0f;
    for (int k = 1; k <= npts; k++) {
        float h_k = 2.0f * (k - 1.0f) / (npts - 1.0f) - 1.0f;
        float theta_k = acosf(h_k);
        if (k == 1 || k == npts) {
            phi_k = 0.0f;
        } else {
            phi_k = fmod((phi_k + 3.6f / sqrtf(npts * (1.0f - h_k*h_k))), (float)(2.0f*M_PI));
        }
        h_points[k-1].x = cosf(phi_k) * sinf(theta_k);
        h_points[k-1].y = sinf(phi_k) * sinf(theta_k);
        h_points[k-1].z = cosf(theta_k);
    }

    if (d_pos == NULL || d_radius == NULL || d_sasa == NULL || d_points == NULL || h_points == NULL || h_sasa == NULL) {
        std::cout << "\033[31m***" << std::endl
            << "*** Error - Allocation of Memory failed!!!" << std::endl
            << "***\033[0m" << std::endl;
    }

    cudaMemcpy(d_pos, pos, static_cast<size_t>(natoms * sizeof(*d_pos)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_radius, radius, static_cast<size_t>(natoms * sizeof(*d_radius)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sasa, h_sasa, static_cast<size_t>(natoms * sizeof(*d_sasa)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_points, h_points, static_cast<size_t>(npts * sizeof(*d_points)), cudaMemcpyHostToDevice);

    SasaKernel<<< gridSize, blockSize >>>(natoms, pairdist, npts, srad, d_pos, d_radius, d_points, d_sasa);
    cudaDeviceSynchronize();

    cudaError_t cudaError = cudaGetLastError();
    if ( cudaError != cudaSuccess ) {
        std::cout << "\033[31m***" << std::endl
                  << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                    << std::endl
                  << "***\033[0m" << std::endl;
    }

    cudaMemcpy(h_sasa, d_sasa, static_cast<size_t>(natoms * sizeof(*d_sasa)), cudaMemcpyDeviceToHost);
    float sasa = h_sasa[5];
    free(h_sasa);
    free(h_points);
    cudaFree(d_sasa);
    cudaFree(d_pos);
    cudaFree(d_radius);
    cudaFree(d_points);


    return sasa;
}
