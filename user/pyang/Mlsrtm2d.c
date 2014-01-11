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
static float *rr, *mm, *gm, *gr, *sm;
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
	mm=(float*)malloc(nm*sizeof(float));
	gm=(float*)malloc(nm*sizeof(float));
	sm=(float*)malloc(nm*sizeof(float));
}

void lsrtm2d_close()
{
	free(rr);
	free(gr);
	free(mm);
	free(gm);
	free(sm);
}

void lsrtm2d(float dz, float dx, float dt, int n0, int n1, 
int n2, int nb, int nt, float **vv, float *mod, float *dat, int niter)
/*< LSRTM with conjugate gradient method >*/
{
	float beta, alpha, g0, gn, gnp, res;

	for(int i=0; i<nd;i++) rr[i]=-dat[i];
	memset(gr,0,nd*sizeof(float));
	memset(mm,0,nm*sizeof(float));
	memset(gm,0,nm*sizeof(float));
	memset(sm,0,nm*sizeof(float));
	rtm2d_init(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat);

	for(int iter=0;iter<niter;iter++)
	{
		rtm2d_lop(true, false, nm, nd, gm, rr); 
		gn=cblas_dsdot(nm,gm,1,gm,1);
		if (iter==0){
			beta=0.0;
			g0=gn;
		}else{
			beta=gn/gnp;
			if(beta<tol || gn/g0<tol) break;
		}
		gnp=gn;

		for(int i=0; i<nm; i++)	sm[i]=gm[i]+beta*sm[i];
		rtm2d_lop(false, false, nm, nd, sm, gr);

		alpha=-gn/cblas_dsdot(nd,gr,1,gr,1);

		for(int i=0; i<nm; i++) mm[i]+=alpha*sm[i];
		for(int i=0; i<nd; i++) rr[i]+=alpha*gr[i];

		res=cblas_dsdot(nd,rr,1,rr,1);
		if (verb) sf_warning("iteration %d; res %g",iter+1, res);
	}
	for(int i=0; i<nm; i++)	mod[i]=mm[i];

	rtm2d_close();
}

int main(int argc, char* argv[])
{   
	bool verb;
    	int niter, n1, n2, nb, nt,n0, nx;
    	float tol, dt,dx,dz, o1,o2;
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
    	if (!sf_getint("nb",&nb)) nb=50;
	/* number (thickness) of ABC boundary grid on each side */
    	if (!sf_getint("n0",&n0)) n0=nb;
	/* shot depth in the grid */
    	if (!sf_getbool("verb",&verb)) verb=false;
	/* verbosity */
    	if (!sf_getint("niter",&niter)) niter=10;
	/* totol number of least-squares iteration*/
    	if (!sf_getfloat("tol",&tol)) tol=1.e-6;
	/* tolerance of inversion */
   
	if (!sf_histint(data,"n1",&nt)) sf_error("n1");
	/* number of time steps */
	if (!sf_histfloat(data,"d1",&dt)) sf_error("d1");
	/* time sampling interval: dt */
	if (!sf_histint(data,"n2",&nx) || nx != n2) 
	sf_error("Need n2=%d in data",n2);
	sf_putint(imag,"n1",n1);
	sf_putint(imag,"n2",n2);
	sf_putfloat(imag,"d1",dz);
	sf_putfloat(imag,"d2",dx);
	sf_putfloat(imag,"o1",o1);
	sf_putfloat(imag,"o2",o2);
	sf_putstring(imag,"label1","Depth");
	sf_putstring(imag,"label2","Distance");

	/* In rtm, vv is the velocity model [modl], which is input parameter; 
	   mod is the image/reflectivity [imag]; dat is seismogram [data]! */
    	vv = sf_floatalloc2(n1,n2);
    	mod = sf_floatalloc(n1*n2);
    	dat = sf_floatalloc(nt*n2);

    	sf_floatread(vv[0],n1*n2,modl);
	sf_floatread(dat,nt*n2,data);

	lsrtm2d_init(n1*n2, nt*n2, tol, verb);
	lsrtm2d(dz, dx, dt, n0, n1, n2, nb, nt, vv, mod, dat, niter);
	lsrtm2d_close();

    	sf_floatwrite(mod,n1*n2,imag);  /* output image */

	free(*vv); 	free(vv);
	free(mod);
	free(dat); 

    	exit(0);
}

