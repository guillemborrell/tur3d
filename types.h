#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>

#ifndef _types_h_
#define _types_h_

typedef struct grid_t{
    int nz, ny, nx;
} grid_t;

typedef struct shape_t{
    double z, y, x;
} shape_t;

int idx(int, int, int, grid_t);

void fourier_coefficients(double *, double *, double *, grid_t, shape_t);

#endif
