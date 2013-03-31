#define NX 32
#define NY 32
#define NZ 32

#include <gsl/gsl_rng.h>
#include "types.h"
#include "io.h"

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

    /* Allocate one of the arrays in physical space */
    double *vort;
    double complex *vort_hat;
    double complex *workzyx;
    double complex *workxyz;

    vort = fftw_alloc_real(phys.nz * phys.ny * phys.nx);
    vort_hat = fftw_alloc_complex(fou.nz * fou.ny * fou.nx);
    workzyx = fftw_alloc_complex(aliased.nz * aliased.ny * aliased.nx);
    workxyz = fftw_alloc_complex(aliased.nz * fou.nx * fou.ny);

    /* Start the program */
    h5err = H5open();

    /* Random number generation for the first field */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    for (i = 0; i< (phys.nz*phys.ny*phys.nx); i++){
        vort[i] = gsl_rng_uniform(r)-0.5;
    }

    gsl_rng_free(r);
    /* End setup test field */

    h5err = dump_file("original.h5", vort, phys);

    /* Plan the 2d fft */
    fftw_plan fft_2d_forward_plan, fft_2d_backward_plan;

    fft_2d_forward_plan = fftw_plan_dft_r2c_2d(phys.ny,phys.nx,vort,workzyx,FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fft_2d_backward_plan = fftw_plan_dft_c2r_2d(phys.ny,phys.nx,workzyx,vort,FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    for(k=0; k<phys.nz;k++){
        fftw_execute_dft_r2c(fft_2d_forward_plan,vort+(k*phys.ny*phys.nx),workzyx+(k*aliased.ny*aliased.nx));
    }
//    fftw_execute_dft_r2c(fft_2d_forward_plan,vort,workzyx);

    /* Filter in fourier to see if everything is aligned as expected */

    for(k=0; k<aliased.nz; k++){
        for(j=0; j<aliased.ny; j++){
            for(i=0; i<aliased.nx; i++){
                if (j > 4 && j < aliased.ny-4){
                    workzyx[idx(k,j,i,aliased)] = 0.0 + 0.0*I;
                }
                if (i > 4){
                    workzyx[idx(k,j,i,aliased)] = 0.0 + 0.0*I;
                }
            }
        }
    }

    for(i=0; i<2*aliased.nx; i++){
        printf("%f, %f\n",creal(workzyx[idx(0,5,i,aliased)]),cimag(workzyx[idx(0,5,i,aliased)]));
    }

    /* Transform back to physical space */
    for(k=0; k<phys.nz;k++){
        fftw_execute_dft_c2r(fft_2d_backward_plan,workzyx+(k*aliased.ny*aliased.nx),vort+(k*phys.ny*phys.nx));
    }
//    fftw_execute_dft_c2r(fft_2d_backward_plan,workzyx,vort);

    h5err = dump_file("result.h5", vort, phys);

    printf("Program done\n");
    fftw_free(vort);
    fftw_free(vort_hat);
    fftw_free(workxyz);
    fftw_free(workzyx);

    h5err = H5close();
    return 0;
}

