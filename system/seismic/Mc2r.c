/* Cartesian-Coordinates to Riemannian-Coordinates interpolation */
/*
  Copyright (C) 2006 Colorado School of Mines
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
#include "c2r.h"

#define LOOPCC(a) for(iz=0;iz<nz;iz++){ for(ix=0;ix<nx;ix++){  {a} }}
#define LOOPRC(a) for(it=0;it<nt;it++){ for(ig=0;ig<ng;ig++){  {a} }}

int main(int argc, char* argv[])
{
    sf_file Fi=NULL, Fo=NULL, Fr=NULL; /* I/O files */
    bool adj, verb, linear;

    sf_axis ax,az,at,ag;
    int ix,iz,it,ig;
    int nx,nz,nt,ng;
    float dx,dz,x0,z0;

    float  **mapCC=NULL,  **mapRC=NULL;
    float ***comCC=NULL, ***comRC=NULL;
    sf_complex **rays=NULL;

    int nn,ii;
    bool comp; /* complex input */

    /* init RSF */
    sf_init(argc,argv);

    Fi = sf_input (  "in");
    Fr = sf_input ("rays");
    Fo = sf_output( "out");

    if(! sf_getbool(  "verb",&verb   ))   verb=false;
    if(! sf_getbool(   "adj",&adj    ))    adj=false;
    if(! sf_getbool("linear",&linear )) linear=true;  

    ag=sf_iaxa(Fr,1); ng=sf_n(ag); if(verb) sf_raxa(ag);
    at=sf_iaxa(Fr,2); nt=sf_n(at); if(verb) sf_raxa(at);

    if(adj) {
	if(! sf_getint  ("a2n",&nz)) nz=1;
	if(! sf_getfloat("a2o",&z0)) z0=0.;
	if(! sf_getfloat("a2d",&dz)) dz=1.;
	az = sf_maxa(nz,z0,dz); sf_setlabel(az,"a2"); if(verb) sf_raxa(az);

	if(! sf_getint  ("a1n",&nx)) nx=1;
	if(! sf_getfloat("a1o",&x0)) x0=0.;
	if(! sf_getfloat("a1d",&dx)) dx=1.;
	ax = sf_maxa(nx,x0,dx); sf_setlabel(ax,"a1"); if(verb) sf_raxa(ax);

	sf_oaxa(Fo,ax,1);
	sf_oaxa(Fo,az,2);
    } else {
	ax = sf_iaxa(Fi,1); nx=sf_n(ax); if(verb) sf_raxa(ax);
	az = sf_iaxa(Fi,2); nz=sf_n(az); if(verb) sf_raxa(az);

	sf_oaxa(Fo,ag,1);
	sf_oaxa(Fo,at,2);
    }

    nn = sf_leftsize(Fi,2);

    rays =sf_complexalloc2(ng,nt);
    sf_complexread(rays[0],ng*nt,Fr);

    c2r_init(ax,az,ag,at,verb);

/*------------------------------------------------------------*/

    comp=false;
    if(SF_COMPLEX == sf_gettype(Fi)) comp=true;

/*------------------------------------------------------------*/

    mapCC=sf_floatalloc2(nx,nz);
    mapRC=sf_floatalloc2(ng,nt);

    if(comp) {
	comCC=sf_floatalloc3(2,nx,nz);
	comRC=sf_floatalloc3(2,ng,nt);

	for(ii=0;ii<nn;ii++) {
	    sf_warning("%d of %d;",ii,nn);
	    if(adj) {
		sf_floatread (comRC[0][0],2*ng*nt,Fi);

		/* REAL */
		LOOPRC( mapRC[it][ig] = comRC[it][ig][0]; );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPCC( comCC[iz][ix][0] = mapCC[iz][ix]; );

		/* IMAGINARY */
		LOOPRC( mapRC[it][ig] = comRC[it][ig][1]; );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPCC( comCC[iz][ix][1] = mapCC[iz][ix]; );

		sf_floatwrite(comCC[0][0],2*nx*nz,Fo);
	    } else {
		sf_floatread (comCC[0][0],2*nx*nz,Fi);
		
		/* REAL */
		LOOPCC( mapCC[iz][ix] = comCC[iz][ix][0]; );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPRC( comRC[it][ig][0] = mapRC[it][ig]; );

		/* IMAGINARY */
		LOOPCC( mapCC[iz][ix] = comCC[iz][ix][1]; );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPRC( comRC[it][ig][1] = mapRC[it][ig]; );

		sf_floatwrite(comRC[0][0],2*ng*nt,Fo);
	    }
	}
	sf_warning(".");

    } else {

	for(ii=0;ii<nn;ii++) {
	    sf_warning("%d of %d;",ii,nn);
	    if(adj) {
		sf_floatread (mapRC[0],ng*nt,Fi);
		c2r(linear,adj,mapCC,mapRC,rays);
		sf_floatwrite(mapCC[0],nx*nz,Fo);
	    } else {
		sf_floatread (mapCC[0],nx*nz,Fi);		
		c2r(linear,adj,mapCC,mapRC,rays);
		sf_floatwrite(mapRC[0],ng*nt,Fo);
	    }
	}
	sf_warning(".");

    }

    exit(0);
}
