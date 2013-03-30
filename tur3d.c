#define NX 32
#define NY 32
#define NZ 32

#include <stdio.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>
#include <hdf5.h>
#include <hdf5_hl.h>


const int DEBUG = 1;

typedef struct grid_t{
    int nz, ny, nx;
} grid_t;

typedef struct shape_t{
    double z, y, x;
} shape_t;

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

void dealias(fftw_complex *in, fftw_complex *out, grid_t aliased, grid_t fou){
    int i,j,k;
    for(k=0; k<fou.nz; k++){
        for(j=0; j<fou.ny; j++){
            for(i=0; i<fou.nx; i++){
                out[idx(k,j,i,fou)] = in[idx(k,j,i,aliased)];
            }
        }
    }
}

void stretch(fftw_complex *in, fftw_complex *out, grid_t fou, grid_t aliased){
    int i,j,k;
    for(k=0; k<fou.nz; k++){
        for(j=0; j<fou.ny; j++){
            for(i=0; i<fou.nx; i++){
                out[idx(k,j,i,aliased)] = in[idx(k,j,i,fou)];
            }
        }
    }
}

int main(int argc, char *argv[]){
    /* hdf5 error code */
    herr_t h5err;

    /* Useful indices */
    int i,j,k;

    /* NX is the base length, NX_phys takes aliasing into account */
    grid_t phys = {NZ*3/2, NY*3/2, NX*3/2};

    /* Dealiased Fourier coefficients*/
    grid_t fou = {NZ/2+1, NY, NX};

    /* Size for aliased Fourier coefficients, maybe not needed */
    grid_t aliased = {phys.nz/2+1, phys.ny, phys.nx};

    /* Geometrical properties of the box */
    shape_t shape = {1.0*M_PI, 1.0*M_PI, 1.0*M_PI};
    shape_t dx = {shape.z/(NZ-1), shape.y/(NY-1), shape.x/(NX-1)};

    /* Create the vectors for the Fourier coefficients. FFTW allocators are used for convenience */
    double *kz, *ky, *kx;
    fourier_coefficients(kz, ky, kx, fou, shape);

    /* Allocate one of the arrays in physical space */
    double *vorticity_r;
    fftw_complex *vorticity_f;
    fftw_complex *vorticity_a;

    vorticity_r = fftw_alloc_real(phys.nz * phys.ny * phys.nx);
    vorticity_a = fftw_alloc_complex(aliased.nz * aliased.ny * aliased.nx);
    vorticity_f = fftw_alloc_complex(fou.nz * fou.ny * fou.nx);

    /* Start the program */
    h5err = H5open();

    /* Create the test vorticity field */
    for(k=0; k<phys.nz; k++){
        for(j=0; j<phys.ny; j++){
            for(i=0; i<phys.nx; i++){
                vorticity_r[idx(k,j,i,phys)] = (double) (k+j);
            }
        }
    }

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    plan_forward = fftw_plan_dft_r2c_3d(phys.nz,phys.ny,phys.nz,
                                        vorticity_r,vorticity_a,
                                        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    plan_backward = fftw_plan_dft_c2r_3d(phys.nz,phys.ny,phys.nz,
                                         vorticity_a,vorticity_r,
                                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    /* Go from phys to fou */
    fftw_execute_dft_r2c(plan_forward,vorticity_r,vorticity_a);

    dealias(vorticity_a,vorticity_f,aliased,fou);
    /* Differentiate in x */
    for(k=0; k<fou.nz; k++){
        for(j=0; j<fou.ny; j++){
            for(i=0; i<fou.nx; i++){
                vorticity_f[idx(k,j,i,fou)] = kx[i]*I*vorticity_f[idx(k,j,i,fou)];
            }
        }
    }
    stretch(vorticity_f,vorticity_a,fou,aliased);

    /* Go from fou to phys */
    fftw_execute_dft_c2r(plan_backward,vorticity_a,vorticity_r);


    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    printf("Program done\n");
    fftw_free(vorticity_r);
    fftw_free(vorticity_f);

    h5err = H5close();
	return 0;
}
