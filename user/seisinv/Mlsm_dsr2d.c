/* 2-D prestack least-squares migration with split-step DSR. */
/*
  Copyright (C) 2012 China University of Petroleum (East China)
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "lsm_dsr2d.h"

int main (int argc, char *argv[])
{
    int nz;		/* depth samples */
    int ny;             /* lateral samples */
    int nx;		/* number of midpoints 	*/
    int nh;             /* number of offsets */
    int nw;		/* number of frequencies */	
    int nt, ntx, nth;   /* boundary taper size */
    int nr;             /* number of reference velocities */
    int npad;           /* padding on offset wavenumber */

    float z0, dz;	/* depth origin, sampling interval */
    float w0, dw;	/* frequency origin, sampling interval */
    float x00, dx;	/* midpoint origin, sampling interval	*/
    float y0, dy;       /* spatial origin, sampling interval	*/
    float h0, dh;       /* offset origin, sampling interval */
    float dt;           /* time error */

    float **slow;

    sf_complex *mod,*dat; /*model vector, data vector */
	sf_complex *x0;       /*initial model value */
    int nmod,ndat;
	int niter;            /*number of iterations */
	float *err;           /* residual vector */
	float *xx;            /* the real part of the estimated model vector */

	int i;

    bool verb;            /* verbosity */
    float eps;            /* dip filter constant          */  
    sf_file in, out, vel,error;

    sf_init(argc,argv);
    in = sf_input("in"); /*input data file*/
    out = sf_output("out"); /*output the real part of estimated model*/
    vel = sf_input("slowness");
	error=sf_output("error"); /*output the niter(number) residuals */

	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype(out,SF_FLOAT);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stability parameter */
	if(!sf_getint("niter",&niter)) niter=10;
	/* number of iterations */

    if (!sf_histint  (in,"n1",&nh)) nh = 1;
    if (!sf_histfloat(in,"d1",&dh)) sf_error ("No d1= in input");
    if (!sf_histfloat(in,"o1",&h0)) h0=0.;

	if (!sf_histint  (in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");
    if (!sf_histfloat(in,"o2",&x00)) x00=0.;

	if (!sf_histint  (in,"n3",&nw)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dw)) sf_error ("No d3= in input");
	if (!sf_histfloat(in,"o3",&w0)) sf_error ("No o3= in input");

    if (!sf_histint  (vel,"n3",&nz)) sf_error ("No n3= in slowness");
    if (!sf_histfloat(vel,"d3",&dz)) sf_error ("No d3= in slowness");
    if (!sf_histfloat(vel,"o3",&z0)) z0=0.; 

    if (!sf_histint  (vel,"n1",&ny)) sf_error ("No n1= in slowness");
    if (!sf_histfloat(vel,"d1",&dy)) sf_error ("No d1= in slowness");
    if (!sf_histfloat(vel,"o1",&y0)) y0=0.;

    if (!sf_getint("nt",&nt)) nt = 1;
    /* taper size */
    ntx = SF_MIN(nt,nx-1);
    nth = SF_MIN(nt,nh-1);
    
    if (!sf_getint("nr",&nr)) nr = 1;
    /* maximum number of references */

    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time error */

    if (!sf_getint("npad",&npad)) npad = 0;
    /* padding on offset wavenumber */

	sf_putint  (out,"n3",nz);
	sf_putfloat(out,"d3",dz);
	sf_putfloat(out,"o3",z0);
	sf_putint  (out,"n2",ny);
	sf_putfloat(out,"d2",dy);
	sf_putfloat(out,"o2",y0);
	sf_putint  (out,"n1",nh);
	sf_putfloat(out,"d1",dh);
	sf_putfloat(out,"o1",h0);

    /* from hertz to radian */
    dw *= 2.*SF_PI; 
    w0 *= 2.*SF_PI;

    slow = sf_floatalloc2(ny,nz);
    sf_floatread(slow[0],ny*nz,vel);
    sf_fileclose(vel);

	ndat=nh*nx*nw;
	dat=sf_complexalloc(ndat);
	nmod=nh*ny*nz;
	mod=sf_complexalloc(nmod);

	/* obtain initial model x0 */
    x0=sf_complexalloc(nmod);    
	for(i=0;i<nmod;i++){
		x0[i]=sf_cmplx(0.0,0.0);
	}

	err=sf_floatalloc(niter);
	xx=sf_floatalloc(nmod);

	sf_warning("nh=%d nx=%d nw=%d ndat=%d   nh=%d ny=%d nz=%d nmod=%d",nh,nx,nw,ndat,nh,ny,nz,nmod);
    
	sf_complexread(dat,ndat,in);
	sf_fileclose(in);

    lsm_dsr2_init(nz,dz, 
	       nh,dh,h0, 
	       nx,dx,x00, 
	       ny,dy,y0,
           nw,dw,w0, 
	       ntx,nth,
	       nr,
	       npad,verb,eps,slow,dt);
    
	/* update model for optimal solution */
    sf_csolver(lsm_dsr2_lop,sf_ccgstep,nmod,ndat,mod,dat,niter,"x0",x0,"err",err,"end");
    
	lsm_dsr2_close();

	/*display residual of every iteration*/
	for(i=0;i<niter;i++)
		sf_warning("%d of %d error: %f",i+1,niter,err[i]);

	for(i=0;i<nmod;i++){
		xx[i]=crealf(mod[i]);
	}

	sf_floatwrite(xx,nmod,out);
	sf_fileclose(out);
    sf_floatwrite(err,niter,error);
	sf_fileclose(error);

    exit (0);
}
