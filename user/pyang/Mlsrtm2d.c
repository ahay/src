/*  2-D zero-offset least-squares reverse time migration (LSRTM)
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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
#include <rsf.h>

#include "rtm2d.h"

static int nm, nd;
static float *rr, *mm, *gm, *gr;
static float tol;
static bool verb;

void lsrtm2d_init(int nm_, int nd_, float tol_, bool verb_)
/*< allocate variables and initialize parameters >*/
{
	nm=nm_;
	nd=nd_;
	tol=tol_;
	verb=verb_;
	rr=(float*)malloc(nd*sizeof(float));
	gr=(float*)malloc(nd*sizeof(float));
	mm=(float*)malloc(nd*sizeof(float));
	gm=(float*)malloc(nd*sizeof(float));
}

void lsrtm2d_close()
/*< free the allocated variables >*/
{
	free(rr);
	free(gr);
	free(mm);
	free(gm);
}



void lsrtm2d(float dz, float dx, float dt, int n0, int n1, 
int n2, int nb, int nt, float **vv, float *mod, float *dat, int niter)
/*< LSRTM with conjugate gradient method >*/
{
  bool forget;
  int i, iter;
	float res0, res;
	for(i=0; i<nd;i++) 	rr[i]=-dat[i];
	memset(gr, 0, nd*sizeof(float));
	memset(mm, 0, nm*sizeof(float));
	memset(gm, 0, nm*sizeof(float));
	rtm2d_init(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat);

	res0=cblas_dsdot(nd, rr, 1, rr, 1);
	for(iter=0;iter<niter;iter++)
	{
		rtm2d_lop(true,  false, nm, nd, gm, rr);// gm=Ft[rr]
		rtm2d_lop(false, false, nm, nd, gm, gr);// gr=F [gm]		
	    	forget = (bool) (0 == (iter+1)%10); // restart every 10 iterations
		/* Claerbout's CG: (mm, rr)=cgstep(mm, rr, gm, gr); */	
		sf_cgstep(forget, nm, nd, mm, gm, rr, gr); 	
		res=cblas_dsdot(nd, rr, 1, rr, 1);
		if (verb) sf_warning("iteration %d; res %g",iter+1, res);
		if (res/res0<tol) break;
	}
	for(i=0; i<nm; i++)	mod[i]=mm[i];

	rtm2d_close();
}

int main(int argc, char* argv[])
{   
	bool verb;
    	int niter, n1, n2, nb, nt, n0, nx;
    	float tol, dt, dx, dz, o1, o2;
    	float *mod, *dat, **vv;      

    	sf_file data, imag, modl;/* I/O files */

    	/* initialize Madagascar */
    	sf_init(argc,argv);

	data = sf_input ("in"); /* seismic data */
	imag = sf_output("out");  /* output image */
    	modl = sf_input ("vel");/* velocity model */
    
    	if (!sf_histint(modl,"n1",&n1)) sf_error("n1");
	/* 1st dimension size */
    	if (!sf_histint(modl,"n2",&n2)) sf_error("n2");
	/* 2nd dimension size */
    	if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
	/* d1 */
    	if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");
	/* d2 */
    	if (!sf_histfloat(modl,"o1",&o1)) sf_error("o1");
	/* o1 */
    	if (!sf_histfloat(modl,"o2",&o2)) sf_error("o2");
	/* o2 */
    	if (!sf_getint("nb",&nb)) nb=20;
	/* number (thickness) of ABC boundary grid on each side */
    	if (!sf_getint("n0",&n0)) n0=0;
	/* shot depth in the grid */
    	if (!sf_getbool("verb",&verb)) verb=false;
	/* verbosity */
    	if (!sf_getint("niter",&niter)) niter=10;
	/* totol number of least-squares iteration*/
    	if (!sf_getfloat("tol",&tol)) tol=1.e-12;
	/* tolerance of inversion */
   
	if (!sf_histint(data,"n1",&nt)) sf_error("n1");
	/* number of time steps */
	if (!sf_histfloat(data,"d1",&dt)) sf_error("d1");
	/* time sampling interval: dt */
	if (!sf_histint(data,"n2",&nx) || nx != n2) 
	sf_error("Need n2=%d in data",n2);
	sf_putint(imag,"n1",n1+2*nb);
	sf_putint(imag,"n2",n2+2*nb);
	sf_putfloat(imag,"d1",dz);
	sf_putfloat(imag,"d2",dx);
	sf_putfloat(imag,"o1",o1-nb*dz);
	sf_putfloat(imag,"o2",o2-nb*dx);
	sf_putstring(imag,"label1","Depth");
	sf_putstring(imag,"label2","Distance");

	/* In rtm, vv is the velocity model [modl], which is input parameter; 
	   mod is the image/reflectivity [imag]; dat is seismogram [data]! */
    	vv = sf_floatalloc2(n1,n2);
    	mod = sf_floatalloc((n1+2*nb)*(n2+2*nb));
    	dat = sf_floatalloc(nt*n2);

    	sf_floatread(vv[0],n1*n2,modl);
	memset(mod,0,(n1+2*nb)*(n2+2*nb)*sizeof(float));
	sf_floatread(dat,nt*n2,data);
/*
// method 1: use my own CG solver, no reweighting
	lsrtm2d_init(n1*n2, nt*n2, tol, verb);
	lsrtm2d(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat, niter);
	lsrtm2d_close();
*/


/* method 2: use bigsolver, no reweighting (=method 1) */
	rtm2d_init(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat);
   	sf_solver(rtm2d_lop, sf_cgstep, (n1+2*nb)*(n2+2*nb), nt*n2, mod, dat, niter, "verb", verb, "end");
	rtm2d_close();



/*
// method 3: IRLS with bigsolver reweighting for L0/L1 sparsity-promotion
    	float *w=sf_floatalloc((n1+2*nb)*(n2+2*nb));
    	for (int i=0; i<(n1+2*nb)*(n2+2*nb); i++) w[i]=1.0f;
	rtm2d_init(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat);
    	for (int iter = 0; iter < niter; iter++) {
   		sf_solver(rtm2d_lop, sf_cgstep, (n1+2*nb)*(n2+2*nb), nt*n2, mod, dat, 1, "x0", mod, "mwt", w, "end");

    		for (int i=0; i<(n1+2*nb)*(n2+2*nb); i++) w[i]=fabsf(mod[i]);//L0-constraint for the model
		//for (int i=0; i<(n1+2*nb)*(n2+2*nb); i++) w[i]=sqrtf(fabsf(mod[i]));//L1-constraint for the model

		if(verb) sf_warning("iteration %d;",iter+1);
    	}
	rtm2d_close();
	free(w);
*/

    	sf_floatwrite(mod, (n1+2*nb)*(n2+2*nb), imag);  /* output image */

	free(*vv); free(vv);
	free(mod);
	free(dat); 

    	exit(0);
}

