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
#include "types.h"


int main(int argc, char *argv[]){
    /* hdf5 error code */
    herr_t h5err;

    /* Useful indices */
    int i,j,k;

    /* NX is the base length, NX_phys takes aliasing into account */
    grid_t phys = {NZ*3/2, NY*3/2, NX*3/2};

    /* Dealiased Fourier coefficients*/
    grid_t fou = {NZ, NY, NX/2+1};

    /* Size for aliased Fourier coefficients, maybe not needed */
    grid_t aliased = {phys.nz, phys.ny, phys.nx/2+1};

    /* Geometrical properties of the box */
    shape_t shape = {1.0*M_PI, 1.0*M_PI, 1.0*M_PI};
    shape_t dx = {shape.z/(NZ-1), shape.y/(NY-1), shape.x/(NX-1)};

    /* Create the vectors for the Fourier coefficients. FFTW allocators are used for convenience */
    double *kz, *ky, *kx;
    fourier_coefficients(kz, ky, kx, fou, shape);

    /* Allocate one of the arrays in physical space */
    double *vorticity_r;
    fftw_complex *vorticity_f;
    fftw_complex *workarray;

    vorticity_r = fftw_alloc_real(phys.nz * phys.ny * phys.nx);
    vorticity_f = fftw_alloc_complex(fou.nz * fou.ny * fou.nx);
    workarray = fftw_alloc_complex(aliased.ny * aliased.nx);

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


//    printf("%f, i%f\n",creal(vorticity_f[0]),cimag(vorticity_f[0]));

//    /* Differentiate in x */
//    for(k=0; k<fou.nz; k++){
//        for(j=0; j<fou.ny; j++){
//            for(i=0; i<fou.nx; i++){
//                vorticity_f[idx(k,j,i,fou)] = -(kx[i]*I) * vorticity_f[idx(k,j,i,fou)];
//            }
//        }
//    }
//    stretch(vorticity_f,vorticity_a,fou,aliased);

//    /* Go from fou to phys */
//    fftw_execute_dft_c2r(plan_backward,vorticity_a,vorticity_r);

    h5err = dump_file("dump.h5", vorticity_r, phys);

//    fftw_destroy_plan(plan_forward);
//    fftw_destroy_plan(plan_backward);

    printf("Program done\n");
    fftw_free(vorticity_r);
    fftw_free(vorticity_f);

    h5err = H5close();
	return 0;
}
