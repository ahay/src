/* Cartesian-Coordinates to Ray-Coordinates interpolation */
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

#include <rsf.h>

#define LOOPCC(a) for(ix=0;ix<ax.n;ix++){ for(iz=0;iz<az.n;iz++){ {a} }}
#define LOOPRC(a) for(ig=0;ig<ag.n;ig++){ for(it=0;it<at.n;it++){ {a} }}

int main(int argc, char* argv[])
{
    sf_file Fi, Fo, Fr; /* I/O files */
    bool adj, verb;

    axa ax,az,at,ag;
    int  ix,iz,it,ig;
    int  jx,jz;
    float x, z;
    float xmin,xmax,zmin,zmax;
    float fx,fz;
    float xo,zo;
    float w00,w01,w10,w11,sumw;

    float **mapCC;
    float **mapRC;
    complex float **rays;
 
    /* init RSF */
    sf_init(argc,argv);

    Fi = sf_input("in");
    Fr = sf_input("rays");
    Fo = sf_output("out");

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getbool( "adj",&adj ))  adj=false;

    iaxa(Fr,&at,1);
    iaxa(Fr,&ag,2);

    raxa(at);
    raxa(ag);

    if(adj) {
	if(! sf_getint  ("azn",&az.n)) az.n=1;
	if(! sf_getfloat("ozn",&az.o)) az.o=0.;
	if(! sf_getfloat("dzn",&az.d)) az.d=1.;

	if(! sf_getint  ("axn",&ax.n)) ax.n=1;
	if(! sf_getfloat("oxn",&ax.o)) ax.o=0.;
	if(! sf_getfloat("dxn",&ax.d)) ax.d=1.;

	oaxa(Fo,&az,1);
	oaxa(Fo,&ax,2);
    } else {
	iaxa(Fi,&az,1);
	iaxa(Fi,&ax,2);

	oaxa(Fo,&at,1);
	oaxa(Fo,&ag,2);
    }

    xmin = ax.o;
    xmax = ax.o + (ax.n-1)*ax.d;
    zmin = az.o;
    zmax = az.o + (az.n-1)*az.d;

    mapCC=sf_floatalloc2  (az.n,ax.n);
    mapRC=sf_floatalloc2  (at.n,ag.n);
    rays =sf_complexalloc2(at.n,ag.n);

    sf_complexread(rays[0],at.n*ag.n,Fr);
    if(adj) {
	sf_floatread(mapRC[0],at.n*ag.n,Fi);
    } else {
	sf_floatread(mapCC[0],az.n*ax.n,Fi);
    }

    for(ig=0;ig<ag.n;ig++) {

	 for(it=0;it<at.n;it++) {

	     z = cimagf(rays[ig][it]);
	     x = crealf(rays[ig][it]);

	     z = SF_MIN(z,zmax);
	     z = SF_MAX(z,zmin);
	     x = SF_MIN(x,xmax);
	     x = SF_MAX(x,xmin);

	     iz=(int)((z-az.o)/az.d);
	     ix=(int)((x-ax.o)/ax.d);

	     iz = SF_MIN(iz  ,az.n);
	     iz = SF_MAX(iz  , 0  );
	     jz = SF_MIN(iz+1,az.n-1);

	     ix = SF_MIN(ix  ,ax.n);
	     ix = SF_MAX(ix  , 0  );
	     jx = SF_MIN(ix+1,ax.n-1);

	     zo = az.o+iz*az.d;
	     xo = ax.o+ix*ax.d;

	     fz= (z-zo)/az.d;
	     fx= (x-xo)/ax.d;
		 
	     fz = SF_MAX(0., fz);
	     fz = SF_MIN(1., fz);
	     fx = SF_MAX(0., fx);
	     fx = SF_MIN(1., fx);

	     /*
	       00   01       --> x
	          .         |
	       10   11      v
	                    z
	     */
	     w00=(1.-fx)*(1.-fz);
	     w10=(   fx)*(1.-fz);
	     w01=(1.-fx)*(   fz);
	     w11=(   fx)*(   fz);

	     /* normalize the sum of the weights */
	     sumw = w00+w10+w01+w11;
	     w00 = w00 / sumw;
	     w10 = w10 / sumw;
	     w01 = w01 / sumw;
	     w11 = w11 / sumw;

	     if(adj) {
		 mapCC[ix][iz] = mapCC[ix][iz] + w00 * mapRC[ig][it];
		 mapCC[ix][jz] = mapCC[ix][jz] + w10 * mapRC[ig][it];
		 mapCC[jx][iz] = mapCC[jx][iz] + w01 * mapRC[ig][it];
		 mapCC[jx][jz] = mapCC[jx][jz] + w11 * mapRC[ig][it];
	     } else {
		 mapRC[ig][it] = mapRC[ig][it] + w00 * mapCC[ix][iz];
		 mapRC[ig][it] = mapRC[ig][it] + w10 * mapCC[ix][jz];
		 mapRC[ig][it] = mapRC[ig][it] + w01 * mapCC[jx][iz];
		 mapRC[ig][it] = mapRC[ig][it] + w11 * mapCC[jx][jz];
	     }

	 } /* it */
    } /* ig */

    if(adj) {
	sf_floatwrite(mapCC[0],az.n*ax.n,Fo);
    } else {
	sf_floatwrite(mapRC[0],at.n*ag.n,Fo);
    }

}
