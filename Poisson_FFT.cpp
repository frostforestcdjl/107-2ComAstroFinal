#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>

#define PI 3.141592653589793



int main( int argc, char *argv[] )
{
    //constant
    char method
    float u
    float r
    float v
    
    
    //FFT solver
    
    const ptrdiff_t Nfft = N;
    fftw_plan plan;
    double *rin;
    fftw_complex *cout;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j, k;
    
    MPI_Init(&argc, &argv);
    fftw_mpi_init();
    
    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_3d(Nfft, Nfft, Nfft/2+1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    rin = fftw_alloc_real(2 * alloc_local);
    cout = fftw_alloc_complex(alloc_local);
    
    /* create plan for out-of-place r2c DFT */
    plan = fftw_mpi_plan_dft_r2c_3d(Nfft, Nfft, Nfft, rin, cout, MPI_COMM_WORLD,
                                    FFTW_MEASURE);
    
    /* initialize rin to some function my_func(x,y,z) */
    for (i = 0; i < local_n0; ++i)
    for (j = 0; j < Nfft; ++j)
    for (k = 0; k < Nfft; ++k)
    rin[(i*Nfft + j) * (2*(Nfft/2+1)) + k] = rho(local_0_start+i, j, k);
    
    /* compute transforms as many times as desired */
    fftw_execute(plan);
    
    fftw_destroy_plan(plan);
    
    MPI_Finalize();
    
    
    float dk = L/N
    D = //
    kx //
    ky //
    kz //
    float DC = D[0][0][0];
    
    //gravity
    float gx[cells][cells][cells];
    float gy[cells][cells][cells];
    float gz[cells][cells][cells];
    float dx = L/cells;
    for (int j=0; j<cells; j++) {
        for (int k=0; k<cells; k++) {
            for (int l=0; l<cells; l++){
                gx[j][k][l] = -(u[j+2][k+1][l+1]-u[j][k+1][l+1])/2/dx;
                gy[j][k][l] = -(u[j+1][k+2][l+1]-u[j+1][k][l+1])/2/dx;
                gz[j][k][l] = -(u[j+1][k+1][l+2]-u[j+1][k+1][l])/2/dx;
            }
        }
    }
    
    //
    float a[3][particle];
    if (method == 'NGP'){
        for (int i=0; i<particle; i++){
            for (int j=0; j<cells; j++) {
                for (int k=0; k<cells; k++) {
                    for (int l=0; l<cells; l++){
                        if (abs(r[0][i] - (1.5*dx + j*dx)) < 0.5*dx && abs(r[1][i] - (1.5*dx + k*dx)) < 0.5*dx && abs(r[2][i] - (1.5*dx + l*dx)) < 0.5*dx){
                            a[0][i] += g[0][j][k][l];
                            a[1][i] += g[1][j][k][l];
                            a[2][i] += g[2][j][k][l];
                        }
                    }
                }
            }
        }
    }
    else if (method == 'CIC'){
        for (int i=0; i<particle; i++){
            for (int j=0; j<cells; j++) {
                for (int k=0; k<cells; k++) {
                    for (int l=0; l<cells; l++){
                        if (abs(r[0][i] - (1.5*dx + j*dx)) < dx && abs(r[1][i] - (1.5*dx + k*dx)) < dx && abs(r[2][i] - (1.5*dx + l*dx)) < dx){
                            a[0][i] += g[0][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx);
                            a[1][i] += g[1][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx);
                            a[2][i] += g[2][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx);
                        }
                    }
                }
            }
        }
    }
    else if (method == 'TSC'){
        for (int i=0; i<particle; i++){
            for (int j=0; j<cells; j++) {
                for (int k=0; k<cells; k++) {
                    for (int l=0; l<cells; l++){
                        if (abs(r[0][i] - (1.5*dx + j*dx)) < (1.5*dx) && abs(r[1][i] - (1.5*dx + k*dx)) < (1.5*dx) && abs(r[2][i] - (1.5*dx + l*dx)) < (1.5*dx)){
                            a[0][i] += g[0][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2));
                            a[1][i] += g[1][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2));
                            a[2][i] += g[2][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2))  * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2));
                        }
                    }
                }
            }
        }
    }
    
    return EXIT_SUCCESS;

    
}
