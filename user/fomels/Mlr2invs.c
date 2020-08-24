/* Post-stack 2D LSRTM Two-step Lowrank with preconditioning */
/*
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

#include <math.h>
#include <rsf.h>

#include "invert.h"
#include "chain2dweights.h"
#include "lr2zortm.h"
#include "fft2.h"



int main(int argc, char* argv[])
{
	int nt, nx, nz, nzpad, nxpad, nk, pad = 1, niter,nz2_wf,nx2_wf,nk_wf;
	int nzx, n2, m2, n3,i3;
	float dt, dx, dz;
	float x0;
	//float t0, x0, z0;
	float *data, *modl, *error=NULL, *ww, *ff, **lft, **rht;
	bool isCmplx = false;
    char *errfile;
    sf_file in, out, left, right, fwght, twght, err=NULL;

    sf_init (argc,argv);

    /* Set up I/O */
    in = sf_input("in"); /* data */
    out = sf_output("out"); /* ls image */

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");

	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",0.);
	sf_putstring(out,"label1","Depth");

	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);
	sf_putstring(out,"label2","Distance");

    /* Time/space and Frequency weights initialization */
	    
	/*
     nk (right.rsf) in low-rank decompose comes from transp of velocity model 
     i.e. transp | fft1 | fft3 axis=2 
     But nk in chain does have transp i.e. fft1 | fft3 axis=2 
     different nk(s) is needed 
     The fact that in low rank example you need to transp velocity model first is to make 
     the the function call nk = fft2_init(...) follows its definition - nx comes before nz
    */

    if (NULL != sf_getstring("fweight")) {
    sf_warning("weights are not NULL !");

    nk_wf = fft2_init(false, 1, nz, nx, &nz2_wf, &nx2_wf);


	
	fwght = sf_input("fweight");
   	ff = sf_floatalloc(nk_wf);
	sf_floatread(ff,nk_wf,fwght);
	sf_fileclose(fwght);

	twght = sf_input("tweight");
	ww = sf_floatalloc(nz*nx);
	sf_floatread(ww,nz*nx,twght);
	sf_fileclose(twght);
	   
	chain2dweights_init(nz,nx,nk_wf,nz2_wf,nx2_wf,ww,ff);

    } else {
	fwght = NULL;
	twght = NULL;
    }

    /* Lowrank propagation parameters */
    nk = fft2_init(isCmplx, pad, nx, nz, &nxpad, &nzpad);
    nzx = nz*nx;
    n3 = sf_leftsize(in,2); 

    /* left and right propagator matrix */
    left = sf_input("left"); 
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("No n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lft = sf_floatalloc2(nzx,m2);
    rht = sf_floatalloc2(m2,nk);

    sf_floatread(lft[0],nzx*m2,left);
    sf_floatread(rht[0],m2*nk,right);


    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */

    /* error file setup */
    if (niter > 0) {
	errfile = sf_getstring("err");
	/* output file for error */
	if (NULL != errfile) {
	    err = sf_output(errfile);
	    sf_putint(err,"n1",niter);
	    sf_putfloat(err,"d1",1);
	    sf_putfloat(err,"o1",1);
	    sf_putstring(err,"label1","Iteration Number");
	    sf_putstring(err,"label2","Relative Squared Error");
	    sf_putint(err,"n2",1);
	}
	error = sf_floatalloc(niter);
    }

    data = sf_floatalloc(nt*nx);
    modl = sf_floatalloc(nz*nx);


	lr2zortm_init(nt, nx, nz, nk, nxpad, nzpad, m2, lft, rht);

    for (i3=0; i3 < n3; i3++) {

		sf_floatread (data,nt*nx,in);

		if (NULL != fwght) {
			sf_warning("With preconditioning !\n ");

		    invert(lr2zortm_lop,chain2dweights_lop,
			   niter,nz*nx,nt*nx,modl,NULL,data,error);
		} else {
			sf_warning("No preconditioning !\n ");
		    invert(lr2zortm_lop,NULL,niter,nz*nx,nt*nx,modl,NULL,data,error);
		}

		sf_floatwrite (modl,nz*nx,out);
		if (NULL != err) sf_floatwrite(error,niter,err);
	}
	exit(0);
}
