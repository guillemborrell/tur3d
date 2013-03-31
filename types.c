#include "stdio.h"
#include "types.h"

int DEBUG = 1;

int idx(int k, int j, int i, grid_t shape){
    return k*shape.ny*shape.nx + j*shape.nx + i;
}

void fourier_coefficients(double *kz, double *ky, double *kx, grid_t fou, shape_t shape){
    int i,j,k;
    int kfac;

    kx = fftw_alloc_real(fou.nx);
    ky = fftw_alloc_real(fou.ny);
    kz = fftw_alloc_real(fou.nz);

    // kx
    if (DEBUG) {printf("kx coefficients\n");}
    for (i=0; i<fou.nx; i++){
        kx[i] = (double)i*(2.0*M_PI)/shape.x;
        if (DEBUG) {printf("%f ", kx[i]);}
    }
    // ky
    if (DEBUG) {printf("\n ky coefficients\n");}
    kfac = fou.ny/2-1;
    for (j=0; j<fou.ny; j++){
        ky[j] = (double)((fou.ny/2+j)%(fou.ny)-fou.ny/2) * (2.0*M_PI)/shape.y;
        if (DEBUG) {printf("%f ", ky[j]);}
    }
    // kz
    if (DEBUG) {printf("\n kz coefficients\n");}
    kfac = fou.nz/2-1;
    for (k=0; k<fou.nz; k++){
        kz[k] = (double)((fou.nz/2+k)%(fou.nz)-fou.nz/2) * (2.0*M_PI)/shape.z;
        if (DEBUG) {printf("%f ", kz[k]);}
    }
    if (DEBUG) {printf("\n");}
}
