/* frequency domain Radon transform */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


//	u = 1/sqrt(pi*Nx) Ld  :  L(np,nx) in Fortran format

#include <rsf.h>

typedef struct tag_frt
{
	int nx,np;			// fart firt fhrt
	double *px, *pp;	// fart firt fhrt
	sf_complex **L;		// fart firt fhrt
	sf_complex **LLH;	// firt fhrt
	int *ipiv;			// firt fhrt
	sf_complex **A, *x;	// fhrt
}frt;


void* sf_frt_init(int mtd, int curv,
	float ox, float dx, int nx,
	float op, float dp, int np)
/*< frt initialize >*/
{
	frt *p;
	int ix, ip;
	double x1;

	p = (frt*)malloc(sizeof(frt));

	p->np = np;
	p->nx = nx;

	p->px = (double*)sf_alloc(nx, sizeof(double));
	p->pp = (double*)sf_alloc(np, sizeof(double));
	p->L  = sf_complexalloc2(np, nx);

	for(ix=0; ix<nx; ix++)
	{
		x1=ox+ix*dx;
		switch(curv)
		{
		case 1:
			p->px[ix] = x1*x1;
			break;
		default:
			p->px[ix] = x1;
		}
	}
	for(ip=0; ip<np; ip++) p->pp[ip] = op+ip*dp;

	if(mtd >= 1)
	{
		p->LLH  = sf_complexalloc2(np, np);
		p->ipiv = sf_intalloc(np);
	}else{
		p->LLH  = NULL;
		p->ipiv = NULL;
	}

	if(mtd >= 2)
	{
		p->A = sf_complexalloc2(np, np);
		p->x = sf_complexalloc(np);
	}else{
		p->A = NULL;
		p->x = NULL;
	}

	return p;

}

void sf_frt_close(void *h)
/*< frt release memory >*/
{
	frt *p = (frt*)h;

	free(p->px);
	free(p->pp);
	free(p->L[0]);
	free(p->L);

	if(p->LLH != NULL)
	{
		free(p->LLH[0]);
		free(p->LLH);
		free(p->ipiv);
		if(p->A != NULL)
		{
			free(p->A[0]);
			free(p->A);
			free(p->x);
		}
	}
	free(p);
}



static void frt_construct(double freq, sf_complex**L, 
		double *px, double *pp, int nx, int np)
{
	double x1, x2, x3;
	int ix, ip;

	x1 = 1.0/sqrt(SF_PI*nx);
	for(ix=0; ix<nx; ix++)
	{
		x2 = 2*SF_PI*freq*px[ix];
		for(ip=0; ip<np; ip++)
		{
			x3 = x2*pp[ip];
			L[ix][ip] = x1*sf_cmplx(cos(x3),sin(x3));
		}
	}


}

void sf_fart(void*h, double freq, sf_complex **in, sf_complex **out, int nn)
/*< forward Radon transform by adjoint operator >*/
{
	sf_complex alpha, beta;
	frt *p;

	p = (frt*) h;
	frt_construct(freq, p->L, p->px, p->pp, p->nx, p->np);

	alpha = 1.0;
	beta  = 0.0;

	cgemm_("N", "N", &p->np, &nn, &p->nx, 
		&alpha, p->L[0], &p->np, in[0], &p->nx,
		&beta, out[0], &p->np);
}

void sf_ifart(void*h, double freq, sf_complex **in, sf_complex **out, int nn)
/*< inverse Radon transform>*/
{
	sf_complex alpha, beta;
	frt *p;

	p = (frt*) h;

	frt_construct(freq, p->L, p->px, p->pp, p->nx, p->np);

	alpha = 1.0;
	beta  = 0.0;
	cgemm_("C", "N", &p->nx, &nn, &p->np, 
		&alpha, p->L[0], &p->np, in[0], &p->np,
		&beta, out[0], &p->nx);
}

void sf_firt(void *h, double freq, 
	sf_complex **in, sf_complex **out, int nn, float mu)
/*< forward Radon transform by LS inversion>*/
{
	int i1, lwork, info;
	sf_complex alpha, beta, *work, workopt;
	frt *p;

	p = (frt*) h;
	frt_construct(freq, p->L, p->px, p->pp, p->nx, p->np);
	

	alpha = 1.0;
	beta  = 0.0;

	cherk_("U", "N", &p->np, &p->nx, (float*) &alpha,
		p->L[0], &p->np, (float*) &beta, p->LLH[0], &p->np);
	for(i1=0; i1<p->np; i1++) p->LLH[i1][i1] += mu;

	cgemm_("N", "N", &p->np, &nn, &p->nx, 
		&alpha, p->L[0], &p->np, in[0], &p->nx,
		&beta, out[0], &p->np);

	lwork = -1;
#if defined(HAVE_MKL)
	chesv_("U", &p->np, &nn, (MKL_Complex8 *)(p->LLH[0]), &p->np, p->ipiv, 
		(MKL_Complex8 *)out[0], &p->np, (MKL_Complex8 *)&workopt, &lwork, &info );
#elif defined(__APPLE__)
	chesv_("U", &p->np, &nn, (__CLPK_complex *)(p->LLH[0]), &p->np, p->ipiv, 
		(__CLPK_complex *)out[0], &p->np, (__CLPK_complex *)&workopt, &lwork, &info );
#else
	chesv_("U", &p->np, &nn, p->LLH[0], &p->np, p->ipiv, 
		out[0], &p->np, &workopt, &lwork, &info );
#endif
	lwork = (int)creal(workopt);
	work = sf_complexalloc(lwork);
#if defined(HAVE_MKL)
	chesv_("U", &p->np, &nn, (MKL_Complex8 *)(p->LLH[0]), &p->np, p->ipiv, 
		(MKL_Complex8 *)out[0], &p->np, (MKL_Complex8 *)work, &lwork, &info );
#elif defined(__APPLE__)
	chesv_("U", &p->np, &nn, (__CLPK_complex *)(p->LLH[0]), &p->np, p->ipiv, 
		(__CLPK_complex *)out[0], &p->np, (__CLPK_complex *)work, &lwork, &info );
#else
	chesv_("U", &p->np, &nn, p->LLH[0], &p->np, p->ipiv, 
		out[0], &p->np, work, &lwork, &info );
#endif
	if(info<0) sf_warning("The %d-th argument had an illegal value", -info);
	if(info>0) sf_warning("Dialogal matrix is singular: %d", info);

	free(work);
}


void sf_fhrt_reg(sf_complex **io, int n1, int n2, float mu, float eta)
/*< caculate regularization from solution >*/
{
	int i1, i2;
	double max;
	sf_complex *p;

	for(i2=0; i2<n2; i2++)
	{
		p = io[i2];
		max = 0.0;
		for(i1=0; i1<n1; i1++) 
		{
			p[i1] *= conjf(p[i1]);
			if(crealf(p[i1]) > max) max = crealf(p[i1]);
		}
		if(max > 0.0)
		{
			for(i1=0; i1<n1; i1++) 
				p[i1] = mu/(eta + p[i1]/max);
		}else 
			for(i1=0; i1<n1; i1++) 
				p[i1] = mu;
	}
}

void sf_fhrt(void *h, double freq, 
	sf_complex **in, sf_complex **out, int nn, 
	float mu, float eta, int niter)
/*< forward high-resolution Radon transform >*/
{
	int lwork, info, nrhs;
	int i1, i2, it;
	sf_complex alpha, beta, *work, workopt;
	frt *p;

	p = (frt*) h;
	frt_construct(freq, p->L, p->px, p->pp, p->nx, p->np);

	alpha = 1.0;
	beta  = 0.0;
	nrhs  = 1;

	cherk_("U", "N", &p->np, &p->nx, (float*) &alpha,
		p->L[0], &p->np, (float*)&beta, p->LLH[0], &p->np);

	cgemm_("N", "N", &p->np, &nn, &p->nx, 
		&alpha, p->L[0], &p->np, in[0], &p->nx,
		&beta, out[0], &p->np);

	lwork = -1;
#if defined(HAVE_MKL)
	chesv_("U", &p->np, &nrhs, (MKL_Complex8 *)(p->LLH[0]), &p->np, p->ipiv, 
		(MKL_Complex8 *)out[0], &p->np, (MKL_Complex8 *)&workopt, &lwork, &info );
#elif defined(__APPLE__)
	chesv_("U", &p->np, &nrhs, (__CLPK_complex *)(p->LLH[0]), &p->np, p->ipiv, 
		(__CLPK_complex *)out[0], &p->np, (__CLPK_complex *)&workopt, &lwork, &info );
#else
	chesv_("U", &p->np, &nrhs, p->LLH[0], &p->np, p->ipiv, 
		p->x, &p->np, &workopt, &lwork, &info );
#endif
	lwork = (int)(creal(workopt)+0.5);
	work = sf_complexalloc(lwork);

	for(i2=0; i2<nn; i2++)
	{
		memcpy(p->x, out[i2], p->np*sizeof(sf_complex));
		for(it=0; it<niter; it++)
		{
			sf_fhrt_reg(& p->x, p->np, 1, mu, eta);
			memcpy(p->A[0], p->LLH[0], p->np*p->np*sizeof(sf_complex));
			for(i1=0; i1<p->np; i1++) 
				p->A[i1][i1] += crealf(p->x[i1]);
			memcpy(p->x, out[i2], p->np*sizeof(sf_complex));
#if defined(HAVE_MKL)
			chesv_("U", &p->np, &nrhs, (MKL_Complex8 *)(p->A[0]), &p->np, p->ipiv,
				   (MKL_Complex8 *)(p->x), &p->np, (MKL_Complex8 *)work, &lwork, &info);
#elif defined(__APPLE__)
			chesv_("U", &p->np, &nrhs, (__CLPK_complex *)(p->A[0]), &p->np, p->ipiv,
				   (__CLPK_complex *)(p->x), &p->np, (__CLPK_complex *)work, &lwork, &info);
#else
			chesv_("U", &p->np, &nrhs, p->A[0], &p->np, p->ipiv,
				   p->x, &p->np, work, &lwork, &info);
#endif
		if (info < 0)
			sf_warning("The %d-th argument had an illegal value", -info);
		if (info > 0)
			sf_warning("Dialogal matrix is singular: %d", info);
		}
		memcpy(out[i2], p->x, p->np*sizeof(sf_complex));
	}
	free(work);
}


void sf_fcrt(void*h, double freq,
	sf_complex **in, sf_complex **out, sf_complex **ref, int nn)
/*< forward frequency cooperative high-resolution Radon transform >*/
{
	int lwork, info, nrhs;
	int i1, i2;
	sf_complex alpha, beta, *work, workopt;
	frt *p;

	p = (frt*) h;
	frt_construct(freq, p->L, p->px, p->pp, p->nx, p->np);

	alpha = 1.0;
	beta  = 0.0;
	nrhs  = 1;

	cherk_("U", "N", &p->np, &p->nx, (float*)&alpha,
		p->L[0], &p->np, (float*)&beta, p->LLH[0], &p->np);

	cgemm_("N", "N", &p->np, &nn, &p->nx, 
		&alpha, p->L[0], &p->np, in[0], &p->nx,
		&beta, out[0], &p->np);

	lwork = -1;
#if defined(HAVE_MKL)
	chesv_("U", &p->np, &nrhs, (MKL_Complex8 *)(p->LLH[0]), &p->np, p->ipiv, 
		(MKL_Complex8 *)(p->x), &p->np, (MKL_Complex8 *)&workopt, &lwork, &info );
#elif defined(__APPLE__)
	chesv_("U", &p->np, &nrhs, (__CLPK_complex *)(p->LLH[0]), &p->np, p->ipiv, 
		(__CLPK_complex *)(p->x), &p->np, (__CLPK_complex *)&workopt, &lwork, &info );
#else
	chesv_("U", &p->np, &nrhs, p->LLH[0], &p->np, p->ipiv, 
		p->x, &p->np, &workopt, &lwork, &info );
#endif
	lwork = (int)(creal(workopt)+0.5);
	work = sf_complexalloc(lwork);

	for(i2=0; i2<nn; i2++)
	{
		memcpy(p->A[0], p->LLH[0], p->np*p->np*sizeof(sf_complex));
		for(i1=0; i1<p->np; i1++) 
			p->A[i1][i1] += crealf(ref[i2][i1]);
#if defined(HAVE_MKL)
		chesv_("U", &p->np, &nrhs, (MKL_Complex8 *)(p->A[0]), &p->np, p->ipiv,
			   (MKL_Complex8 *)out[i2], &p->np, (MKL_Complex8 *)work, &lwork, &info);
#elif defined(__APPLE__)
		chesv_("U", &p->np, &nrhs, (__CLPK_complex *)(p->A[0]), &p->np, p->ipiv,
			   (__CLPK_complex *)out[i2], &p->np, (__CLPK_complex *)work, &lwork, &info);
#else
		chesv_("U", &p->np, &nrhs, p->A[0], &p->np, p->ipiv, 
			out[i2], &p->np, work, &lwork, &info );
#endif
		if(info<0) 
			sf_warning("The %d-th argument had an illegal value", -info);
		if(info>0) 
			sf_warning("Dialogal matrix is singular: %d", info);
	}
}

