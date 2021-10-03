/* Multicomponent data registration analysis with non-stationary model smoothing. */
/*
  Copyright (C) Liu et al. (2021)
  
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
  
  Reference: Liu, X.., X. Chen, M. Bai, and Y. Chen, 2021, Time-lapse image registration by high-resolution time-shift scan, 86, XX-XX, doi: 10.1190/geo2020-0459.1.
*/

#include <string.h>
#include <math.h>
#include <float.h>

#include <rsf.h> 

#include "warpscann.h"

int main(int argc, char* argv[])
{ 
    bool shift, verb;
    int   *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */
    int box[SF_MAX_DIM];
    int i, i1, mode, b, n1, m[4], ntr, n2, order, ng, rect[4], niter, n2g, dim;
    char key[8];
    float **inp, **oth, o1, d1, o2, d2, g0, dg, o, d, *rat1;
    sf_file in, warped, other;
	
    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m[2])) m[2] = 1;
    if(!sf_histint(in,"n3",&m[3])) m[3] = 1;
    ntr = m[2]*m[3];

    if (m[3] > 1) {
	dim = 4;
    } else if (m[2] > 1) {
	dim = 3;
    } else {
	dim = 2;
    }

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getint("ng",&ng)) ng=1;
    /* number of gamma values */
    if (!sf_getfloat("g0",&g0)) sf_error("Need g0=");
    /* gamma origin */
    if (!sf_getfloat("dg",&dg)) dg=g0;
    /* gamma sampling */

    if (!sf_getint("mode",&mode)) mode=0;
    /* mode=0: traditional; mode=1: high-resolution */
    
    if(mode==0)
    {
    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* vertical smoothing */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* gamma smoothing */
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* in-line smoothing */
    if (!sf_getint("rect4",&rect[3])) rect[3]=1;
    /* cross-line smoothing */
    }
    
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */
    if (!sf_getbool("shift",&shift)) shift=false;
    /* use shift instead of stretch */
    
    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;
    sf_putint(warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    m[0] = n2;
    m[1] = ng;

    sf_putint  (warped,"n2",ng);
    sf_putfloat(warped,"d2",dg);
    sf_putfloat(warped,"o2",g0);

    if (dim > 3) {
	if(!sf_histfloat(in,"d3",&d)) d=1.;
	if(!sf_histfloat(in,"o3",&o)) o=0.;
	
	sf_putint  (warped,"n4",m[3]);
	sf_putfloat(warped,"d4",d);
	sf_putfloat(warped,"o4",o);
    } 

    if (dim > 2) {
	if(!sf_histfloat(in,"d2",&d)) d=1.;
	if(!sf_histfloat(in,"o2",&o)) o=0.;
	
	sf_putint  (warped,"n3",m[2]);
	sf_putfloat(warped,"d3",d);
	sf_putfloat(warped,"o3",o);
    }

    if(!sf_getint("accuracy",&order)) {
	/* [1-4] interpolation accuracy */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;
    
    n2g = n2*ng*ntr;
    
    inp = sf_floatalloc2 (n1,ntr);
    oth = sf_floatalloc2 (n2,ntr);
    rat1 = sf_floatalloc (n2g);


    if(mode==1)
    {
	sf_file RCT[SF_MAX_DIM], SFT[SF_MAX_DIM]; 
	/*Calculate dim*/
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    RCT[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(RCT[i])) sf_error("Need float %s",key);
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		SFT[i] = sf_input(key);
		if (SF_INT != sf_gettype(SFT[i])) sf_error("Need int %s",key);
	    } else {
		SFT[i] = NULL;
	    }
	} else {
	    RCT[i] = NULL;
	    SFT[i] = NULL;
	}
    }

	/*reading the non-stationary smoothing radii*/
    for (i=0; i < dim; i++) {
	box[i] = 1;
	if (NULL != RCT[i]) {
	    rct[i] = sf_floatalloc (n2g);
	    sft[i] = sf_intalloc (n2g);
	    sf_floatread(rct[i],n2g,RCT[i]);
	    sf_fileclose(RCT[i]);
	    if (NULL != SFT[i]) {
		sf_intread(sft[i],n2g,SFT[i]);
		sf_fileclose(SFT[i]);
	    } else {
		for (i1=0; i1 < n2g; i1++) {
		    sft[i][i1] = 0;
		}
	    }
	    for (i1=0; i1 < n2g; i1++) {
		b = ceilf(rct[i][i1])+SF_ABS(sft[i][i1]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }  
    }

    if(mode==0)
    {
    	sf_warning("Traditional Time-shift Analysis");
    	warpscan_init(n1,o1,d1,n2,o2,d2,ng,g0,dg,ntr,order,dim,m,rect,niter,shift,verb);
    }
    else
    {
    	sf_warning("High-resolution Time-shift Analysis");
    	sf_warning("dim=%d,m[0]=%d,m[1]=%d,m[2]=%d,m[3]=%d",dim,m[0],m[1],m[2],m[3]);
    	sf_warning("box[0]=%d,box[1]=%d,box[2]=%d,box[3]=%d",box[0],box[1],box[2],box[3]);
    	warpscann_init(n1,o1,d1,n2,o2,d2,ng,g0,dg,ntr,order,dim,m,box,rct,sft,niter,shift,verb);
    }

    sf_floatread(inp[0],n1*ntr,in);
    sf_floatread(oth[0],n2*ntr,other);

    if(mode==0)
    {
    	warpscan(inp,oth,rat1);
	}else{
    	warpscann(inp,oth,rat1);
	}
	
    sf_floatwrite(rat1,n2g,warped);

    exit (0);
}

