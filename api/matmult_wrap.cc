#include "matmult.hh"
#include "matmult_wrap.hh"

#include <rsf.h>

static MatMult oper(0,0,NULL);
// put an object in static storage

void matmult_init (int nx, int ny, float** bb) 
// initialize everything needed
{
    oper = MatMult(nx,ny,bb);
}

void matmult_lop (bool adj, bool add, int nx, int ny, float* x, float*y)
// access it with a C-wrapped interface
{
    if (adj) {
	oper.adjoint(add,x,y);
    } else {
	oper.forward(add,x,y);
    }
}
