/* Conjugate-gradient with shaping regularization using CUDA */
/*
  Copyright (C) 2014 Xi'an Jiaotong University, Pengliang Yang
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "cuda_def.cu"

static int np, nx, nr, nd;
static float *xx, *dd, *r, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

void cuda_conjgrad_init(int np1     /* preconditioned size */, 
		      int nx1     /* model size */, 
		      int nd1     /* data size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */) 
/*< solver constructor >*/
{
    np = np1; 
    nx = nx1;
    nr = nr1;
    nd = nd1;
    eps = eps1*eps1;
    tol = tol1;
    verb = verb1;
    hasp0 = hasp01;

    cudaMalloc(&xx, nx*sizeof(float));
    cudaMalloc(&dd, nd*sizeof(float));
    cudaMalloc(&r, nr*sizeof(float));
    cudaMalloc(&sp, np*sizeof(float));
    cudaMalloc(&gp, np*sizeof(float));
    cudaMalloc(&sx, nx*sizeof(float));
    cudaMalloc(&gx, nx*sizeof(float));
    cudaMalloc(&sr, nr*sizeof(float));
    cudaMalloc(&gr, nr*sizeof(float));
}

void cuda_conjgrad_close(void) 
/*< Free allocated space >*/
{
    cudaFree(xx);
    cudaFree(dd);
    cudaFree(r);
    cudaFree(sp);
    cudaFree(gp);
    cudaFree(sx);
    cudaFree(gx);
    cudaFree(sr);
    cudaFree(gr);
}

void cuda_conjgrad(sf_operator prec  /* data preconditioning */, 
		 sf_operator oper  /* linear operator */, 
		 sf_operator shape /* shaping operator */, 
		 float* p          /* preconditioned model */, 
		 float* x          /* estimated model */, 
		 float* dat        /* data */, 
		 int niter         /* number of iterations */) 
/*< Conjugate gradient solver with shaping >*/
{   
    cudaMemcpy(xx, x, nx*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dd, dat, nd*sizeof(float), cudaMemcpyHostToDevice);

    float gn, gnp, alpha, beta, g0, dg, r0;
    float beta1, beta2, beta3, normr;
    float *d=NULL, *ptr;
    int iter;
    
    if (NULL != prec) {
	cudaMalloc(&d, nd*sizeof(float));
	cuda_scale<<<(nd+Block_Size-1)/Block_Size,Block_Size>>>(d, dat, -1, nd);
	prec(false,false,nd,nr,d,r);
    } else {
	cuda_scale<<<(nr+Block_Size-1)/Block_Size,Block_Size>>>(r, dat, -1, nr);
    }
    
    if (hasp0) { /* initial p */
	shape(false,false,np,nx,p,x);
	if (NULL != prec) {
	    oper(false,false,nx,nd,x,d);
	    prec(false,true,nd,nr,d,r);
	} else {
	    oper(false,true,nx,nr,x,r);
	}
    } else {
	cudaMemset(p, 0, np*sizeof(float));
	cudaMemset(x, 0, nx*sizeof(float));
    } 
    
    dg = g0 = gnp = 0.;
    cuda_dot<<<1,Block_Size>>>(r, r, nr, &r0);

    if (r0 == 0.) {
	if (verb) sf_warning("zero residual: r0=%g",r0);
	return;
    }

    for (iter=0; iter < niter; iter++) {
	cuda_scale<<<(np+Block_Size-1)/Block_Size,Block_Size>>>(gp, p, eps, np);
	cuda_scale<<<(nx+Block_Size-1)/Block_Size,Block_Size>>>(gx, x, -eps, nx);
	if (NULL != prec) {
	    prec(true,false,nd,nr,d,r);
	    oper(true,true,nx,nd,gx,d);
	} else {
	    oper(true,true,nx,nr,gx,r);
	}

	shape(true,true,np,nx,gp,gx);
	shape(false,false,np,nx,gp,gx);

	if (NULL != prec) {
	    oper(false,false,nx,nd,gx,d);
	    prec(false,false,nd,nr,d,gr);
	} else {
	    oper(false,false,nx,nr,gx,gr);
	}

    	cuda_dot<<<1, Block_Size>>>(gp, gp, np, &gn);

	if (iter==0) {
	    g0 = gn;
    	    cudaMemcpy(sp, gp, np*sizeof(float), cudaMemcpyDeviceToDevice);
	    cudaMemcpy(sx, gx, nx*sizeof(float), cudaMemcpyDeviceToDevice);
	    cudaMemcpy(sr, gr, nr*sizeof(float), cudaMemcpyDeviceToDevice);
	} else {
	    alpha = gn / gnp;
	    dg = gn / g0;

	    if (alpha < tol || dg < tol) {
		if (verb) 
		    sf_warning(
			"convergence in %d iterations, alpha=%g, gd=%g",
			iter,alpha,dg);
		break;
	    }

    	    cuda_axpy<<<(np+Block_Size-1)/Block_Size,Block_Size>>>(sp, gp, alpha, np);
	    ptr=sp; sp=gp; gp=ptr;

    	    cuda_axpy<<<(nx+Block_Size-1)/Block_Size,Block_Size>>>(sx, gx, alpha, nx);
	    ptr=sx; sx=gx; gx=ptr;

    	    cuda_axpy<<<(nr+Block_Size-1)/Block_Size,Block_Size>>>(sr, gr, alpha, nr);
	    ptr=sr; sr=gr; gr=ptr;
	}

    	cuda_dot<<<1,Block_Size>>>(sr, sr, nr, &beta1);
    	cuda_dot<<<1,Block_Size>>>(sp, sp, np, &beta2);
    	cuda_dot<<<1,Block_Size>>>(sx, sx, nx, &beta3);

	beta = beta1 + eps*(beta2-beta3);
    	cuda_dot<<<1,Block_Size>>>(r, r, nr, &normr);

	if (verb) sf_warning("iteration %d res: %f grad: %f",
			     iter,normr/r0,dg);

	alpha = - gn / beta;

    	cuda_axpy<<<(np+Block_Size-1)/Block_Size,Block_Size>>>(sp, p, alpha, np);
    	cuda_axpy<<<(nx+Block_Size-1)/Block_Size,Block_Size>>>(sx, x, alpha, nx);
    	cuda_axpy<<<(nr+Block_Size-1)/Block_Size,Block_Size>>>(sr, r, alpha, nr);

	gnp = gn;
    }

    if (NULL != prec) cudaFree(d);

    cudaMemcpy(x, xx, nx*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(dat, dd, nd*sizeof(float), cudaMemcpyDeviceToHost);
}
