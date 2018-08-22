/* Anisotropic diffusion by regularized inversion. Instead of a gradient PWDs in inline and crossline directions are used. 3D. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <assert.h>

#include <rsf.h>

#include "pwddiffuse.h"

int main(int argc, char* argv[])
{
    int niter, i, repeat, i1, i2, i3, n12, n123, n23, n1,n2,n3;
    int nw, nj1, nj2;
    bool sm, adj, test;
    float *data, *data2, eps, *vx, *vy;
    float *pp1, *pp2, *pwddata;
    sf_file in, out, fvx, fvy, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input ("dip");
    fvx = sf_input("vx");
    fvy = sf_input("vy");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(dip) ||
	SF_FLOAT != sf_gettype(fvx) ||
	SF_FLOAT != sf_gettype(fvy)) sf_error("Need float type");

    /* Get dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2=");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3=");

    /* PWD parameters */
    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
	
    if (nw < 1 || nw > 3) sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing iline */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing xline */

    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform PWD filtering */

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of conjugate-gradient iterations */
    if (!sf_getint("repeat",&repeat)) repeat=1;
    /* number of smoothing iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */

    
    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag - when test=y */
    if (!sf_getbool("test",&test)) test=false;
    /* test - applied in either forward or adjoint mode (no inversion) */

    n12 = n1*n2;
    n123 = n12*n3;
    n23 = n2*n3;

    /* Allocate and read data */
    data = sf_floatalloc(n123);
    data2 = sf_floatalloc(n123);
    sf_floatread(data,n123,in);

    /* Allocate space for dip */
    pp1 = sf_floatalloc(n123);
    pp2 = sf_floatalloc(n123);

    /* reading iline dip */
    sf_floatread(pp1,n123,dip);
    
    /* reading xline dip */
    sf_floatread(pp2,n123,dip);

    /* Allocate and read orientation vector */
    vx = sf_floatalloc(n123);
    sf_floatread(vx,n123,fvx);
    vy = sf_floatalloc(n123);
    sf_floatread(vy,n123,fvy);

    pwddiffuse_init(n1,n2,n3 /* data size */,
                    pp1,pp2 /* inline and crossline dips */, 
	            nw,nj1,nj2 /* PWD parameters */,
	            vx,vy /* parallel to edges vector components */,
                    niter /* number of iterations */,
                    repeat /* number of smoothing iterations */,
                    eps /* regularization */);

    if(adj==false && test){
	
	pwddiffuse_fwdoradj(adj,false,n123,n123,data,data2);
    
    }

    if(adj==true && test){

	pwddiffuse_fwdoradj(adj,false,n123,n123,data2,data);

    }

    if (!test) {

	pwddiffuse_lop(n123,n123,data,data2);

    }

    sf_floatwrite(data2,n123,out);

    exit(0);
}
