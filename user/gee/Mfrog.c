/* Simple 2-D wave propagation */
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
#include <stdio.h>

#include <rsf.h>

#include "laplacian.h"

int main(int argc, char* argv[])
{
    bool verb;           
    int type;
    int it,iz,ix;        /* index variables */
    int nt,nz,nx;
    float dt,dz,dx,dt2,old;

    float  *ww,**vv,**rr; /* I/O arrays*/
    float **um,**uo,**ud; /* tmp arrays */

    sf_file Fw,Fv,Fr,Fo; /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(! sf_getint("type", &type)) type=0; /* Laplacian type */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fv = sf_input ("vel");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);

    /* read wavelet, velocity & reflectivity */
    ww=sf_floatalloc(nt);     sf_floatread(ww   ,nt   ,Fw);
    vv=sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fv);
    rr=sf_floatalloc2(nz,nx); sf_floatread(rr[0],nz*nx,Fr);

    /* allocate temporary arrays */
    um=sf_floatalloc2(nz,nx);
    uo=sf_floatalloc2(nz,nx);
    ud=sf_floatalloc2(nz,nx);
 
    dt2 = dt*dt;   
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    ud[ix][iz]=0;
	    vv[ix][iz] *= vv[ix][iz]*dt2;
	}
    }
    laplacian_init(type,nz,nx,dz,dx,vv);

    /* MAIN LOOP */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

	laplacian(uo,ud);

	for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
		/* scale by velocity */
		ud[ix][iz] *= vv[ix][iz];

		/* inject wavelet */
		ud[ix][iz] += ww[it] * rr[ix][iz];
	
		/* time step */
		old = uo[ix][iz];
		uo[ix][iz] += uo[ix][iz] - um[ix][iz] + ud[ix][iz]; 		
		um[ix][iz] = old;
	    }
	}
	
	/* write wavefield to output */
	sf_floatwrite(uo[0],nz*nx,Fo);
    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
