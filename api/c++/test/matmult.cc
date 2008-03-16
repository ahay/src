#include "matmult.hh"

// constructor
MatMult::MatMult(int nx, int ny, float **matrix)
{
    matrix_ = matrix;
    nx_ = nx;
    ny_ = ny;
}

void 
MatMult::forward(bool add, const float *x, float *y)
{
    int ix, iy;

    for (iy = 0; iy < ny_; iy++) {
	if (!add) y[iy] = 0.;
	for (ix = 0; ix < nx_; ix++) {
	    y[iy] += matrix_[iy][ix] * x[ix];
	}
    }
}

void 
MatMult::adjoint(bool add, float *x, const float *y)
{
    int ix, iy;

    for (ix = 0; ix < nx_; ix++) {
	if (!add) x[ix] = 0.;
	for (iy = 0; iy < ny_; iy++) {
	    x[ix] += matrix_[iy][ix] * y[iy];
	}
    }
}
