/* Error by substituting numerical solution into equation */
/*
   Equation for escape positions or angle is (with appropriate BC)
   -S.sin(a).dX/dx - S.cos(a).dX/dz - [cos(a).Sx - sin(a).Sz].dX/da = 0

   Equation for escape traveltime T is (with appropriate BC)
   -S.sin(a).dT/dx - S.cos(a).dT/dz - [cos(a).Sx - sin(a).Sz].dT/da = -S*S
*/
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <assert.h>

#include "gsray.h"

int main(int argc, char* argv[])
{
    int nz,nx,na;
    float oz,ox,oa;
    float dz,dx,da;

    int ix0, ix1, ixs;
    int iz0, iz1, izs;
 
    int iq;                        /* escape variables switch */  

    int ix;                        /* grid points in x */
    int iz;                        /* grid points in z */
    int ia;                        /* grid points in a */

    float a;                       /* angle */

    float ***t;                    /* escape variable */
    float ***err;                  /* substitution error */ 
    float **s,**sx,**sz;           /* slowness, gradients */

    float new_val;

    float ss,ssx,ssz;
    float cs,sn;

    sf_file in,out,slow,slowz,slowx;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    float err_cutoff;

    if (!sf_getfloat("err_cutoff",&err_cutoff)) err_cutoff = 0.2;

    if (!sf_getint("iq",&iq)) iq=2;
    /* switch for escape variable 0=x, 1=a, 2=t, 3=z */

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    /* angle in degrees */
    if (!sf_histint(in,"n3",&na)) sf_error("No n3");
    if (!sf_histfloat(in,"d3",&da)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&oa)) sf_error("No o3=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&oz)) sf_error("No o1=");

    sf_putint (out, "n1", nz);
    sf_putint (out, "n2", nx);
    sf_putint (out, "n3", na);
  
    /* memory allocations */
    t = sf_floatalloc3(nz,nx,na);
    err = sf_floatalloc3(nz,nx,na);

    /* read input escape variable */
    sf_floatread(t[0][0],na*nx*nz,in);

    /* read auxiliary slowness file */
    slow = sf_input("slow");
    s = sf_floatalloc2(nz,nx);
    sf_floatread(s[0],nz*nx,slow);

    /* read auxiliary slowness z-gradient file */
    slowz = sf_input("slowz");
    sz = sf_floatalloc2(nz,nx);
    sf_floatread(sz[0],nz*nx,slowz);

    /* read auxiliary slowness x-gradient file */
    slowx = sf_input("slowx");
    sx = sf_floatalloc2(nz,nx);
    sf_floatread(sx[0],nz*nx,slowx);

    /* convert to radians */
    oa *= SF_PI/180.;
    da *= SF_PI/180.;

    gsray_init(nz,nx,na,
	       oz,ox,oa,
	       dz,dx,da);

    /* memset(err[0][0],(int)0, nz*nx*na*sizeof(float)); */

    for (ia = 0; ia < na; ia++) {

	a = oa + ia*da;

	cs = cosf(a);
	sn = sinf(a);

	if (1e-6 < -sn) {
	    
	    ix0=nx-2; ix1=-1; ixs=-1;
	    for (iz = 0; iz < nz; iz++) 
		boundary_mat(err,iz,nx-1,ia,2);

	} else {
	    
	    ix0=1; ix1=nx; ixs=1;
	    for (iz = 0; iz < nz; iz++) 
		boundary_mat(err,iz,0,ia,2);
	    
	}

	if (1e-6 < -cs) {
	    
	    iz0=nz-2; iz1=-1; izs=-1;
	    for (ix = 0; ix < nx; ix++) 
		boundary_mat(err,nz-1,ix,ia,2);
	    
	} else {
	    
	    iz0=1; iz1=nz; izs=1;
	    for (ix = 0; ix < nx; ix++) 
		boundary_mat(err,0,ix,ia,2);
	    
	}


	for (ix = ix0; ix != ix1; ix += ixs) {
	    
	    for (iz = iz0; iz != iz1; iz += izs) {

		ss = s[ix][iz];
		ssx = sx[ix][iz];
		ssz = sz[ix][iz];

		/* Gauss-Seidel update */			    
		new_val = gs_update(t,cs*ss, sn*ss, (cs*ssx - sn*ssz), iz,ix,ia, ss*ss, ss, 0, iq,
		                    0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

		err[ia][ix][iz] = fabsf(t[ia][ix][iz] - new_val);

	    }/* iz */
	    
	} /* ix */

    } /* ia */

    /* output */
    sf_floatwrite(err[0][0],na*nx*nz,out);
    
    exit(0);
}

