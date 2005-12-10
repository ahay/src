/* Cartesian-Coordinates to Riemannian-Coordinates interpolation */
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

#include "c2r.h"

#define LOOPCC(a) for(iz=0;iz<az.n;iz++){ for(ix=0;ix<ax.n;ix++){  {a} }}
#define LOOPRC(a) for(it=0;it<at.n;it++){ for(ig=0;ig<ag.n;ig++){  {a} }}

int main(int argc, char* argv[])
{
    sf_file Fi, Fo, Fr; /* I/O files */
    bool adj, verb, linear;
    
    axa ax,az,at,ag;
    int ix,iz,it,ig;

    float         **mapCC, **mapRC;
    complex float **comCC, **comRC;
    complex float **rays;

    int nn,ii;
    bool comp; // complex input
 
    /* init RSF */
    sf_init(argc,argv);

    Fi = sf_input (  "in");
    Fr = sf_input ("rays");
    Fo = sf_output( "out");

    if(! sf_getbool(  "verb",&verb   ))   verb=false;
    if(! sf_getbool(   "adj",&adj    ))    adj=false;
    if(! sf_getbool("linear",&linear )) linear=true;  

    iaxa(Fr,&ag,1); if(verb) raxa(ag);
    iaxa(Fr,&at,2); if(verb) raxa(at);

    if(adj) {
	az.l="a2";
	if(! sf_getint  ("a2n",&az.n)) az.n=1;
	if(! sf_getfloat("a2o",&az.o)) az.o=0.;
	if(! sf_getfloat("a2d",&az.d)) az.d=1.;
	if(verb) raxa(az);

	ax.l="a1";
	if(! sf_getint  ("a1n",&ax.n)) ax.n=1;
	if(! sf_getfloat("a1o",&ax.o)) ax.o=0.;
	if(! sf_getfloat("a1d",&ax.d)) ax.d=1.;
	if(verb) raxa(ax);

	oaxa(Fo,&ax,1);
	oaxa(Fo,&az,2);
    } else {
	iaxa(Fi,&ax,1); if(verb) raxa(ax);
	iaxa(Fi,&az,2); if(verb) raxa(az);

	oaxa(Fo,&ag,1);
	oaxa(Fo,&at,2);
    }

    nn = sf_leftsize(Fi,2);

    rays =sf_complexalloc2(ag.n,at.n);
    sf_complexread(rays[0],ag.n*at.n,Fr);
    
    c2r_init(ax,az,ag,at,verb);

/*------------------------------------------------------------*/

    comp=false;
    if(SF_COMPLEX == sf_gettype(Fi)) comp=true;

/*------------------------------------------------------------*/

    mapCC=sf_floatalloc2(ax.n,az.n);
    mapRC=sf_floatalloc2(ag.n,at.n);

    if(comp) {
	comCC=sf_complexalloc2(ax.n,az.n);
	comRC=sf_complexalloc2(ag.n,at.n);

	for(ii=0;ii<nn;ii++) {	
	    if(adj) {
		sf_complexread (comRC[0],ag.n*at.n,Fi);

		// REAL
		LOOPRC( mapRC[it][ig]  = crealf(comRC[it][ig]); );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPCC( comCC[iz][ix]  =        mapCC[iz][ix];  );

		// IMAGINARY
		LOOPRC( mapRC[it][ig]  = cimagf(comRC[it][ig]); );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPCC( comCC[iz][ix] +=      I*mapCC[iz][ix];  );

		sf_complexwrite(comCC[0],ax.n*az.n,Fo);
	    } else {
		sf_complexread (comCC[0],ax.n*az.n,Fi);
		
		// REAL
		LOOPCC( mapCC[iz][ix]  = crealf(comCC[iz][ix]); );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPRC( comRC[it][ig]  =        mapRC[it][ig];  );

		// IMAGINARY
		LOOPCC( mapCC[iz][ix]  = cimagf(comCC[iz][ix]); );
		c2r(linear,adj,mapCC,mapRC,rays);
		LOOPRC( comRC[it][ig] +=      I*mapRC[it][ig];  );

		sf_complexwrite(comRC[0],ag.n*at.n,Fo);
	    }
	}

    } else {	

	for(ii=0;ii<nn;ii++) {	
	    if(adj) {
		sf_floatread (mapRC[0],ag.n*at.n,Fi);
		c2r(linear,adj,mapCC,mapRC,rays);
		sf_floatwrite(mapCC[0],ax.n*az.n,Fo);
	    } else {
		sf_floatread (mapCC[0],ax.n*az.n,Fi);		
		c2r(linear,adj,mapCC,mapRC,rays);
		sf_floatwrite(mapRC[0],ag.n*at.n,Fo);
	    }
	}

    }
    exit(0);
}
