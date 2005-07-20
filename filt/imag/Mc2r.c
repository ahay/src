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

#define LOOPCC(a) for(ix=0;ix<ax.n;ix++){ for(iz=0;iz<az.n;iz++){ {a} }}
#define LOOPRC(a) for(ig=0;ig<ag.n;ig++){ for(it=0;it<at.n;it++){ {a} }}

int main(int argc, char* argv[])
{
    sf_file Fi, Fo, Fr; /* I/O files */
    bool adj, verb, linear;

    axa ax,az,at,ag;
    int  ix,iz,it,ig;
    int  jx,jz;
    float x, z;

    float xmin,xmax,zmin,zmax;
    float fx,fz;
    float xo,zo;
    float w00,w01,w10,w11,sumw;

    float dr,r,sincr;
    int   nsz,nsx;
    float  dz,dx;

    float **mapCC;
    float **mapRC;
    complex float **rays;
 
    /* init RSF */
    sf_init(argc,argv);

    Fi = sf_input("in");
    Fr = sf_input("rays");
    Fo = sf_output("out");

    if(! sf_getbool(  "verb",&verb   ))   verb=false;
    if(! sf_getbool(   "adj",&adj    ))    adj=false;
    if(! sf_getbool("linear",&linear )) linear=true;  

    iaxa(Fr,&at,1); if(verb) raxa(at);
    iaxa(Fr,&ag,2); if(verb) raxa(ag);

    if(adj) {
	az.l="z";
	if(! sf_getint  ("azn",&az.n)) az.n=1;
	if(! sf_getfloat("azo",&az.o)) az.o=0.;
	if(! sf_getfloat("azd",&az.d)) az.d=1.;
	if(verb) raxa(az);

	ax.l="x";
	if(! sf_getint  ("axn",&ax.n)) ax.n=1;
	if(! sf_getfloat("axo",&ax.o)) ax.o=0.;
	if(! sf_getfloat("axd",&ax.d)) ax.d=1.;
	if(verb) raxa(ax);

	oaxa(Fo,&az,1);
	oaxa(Fo,&ax,2);
    } else {
	iaxa(Fi,&az,1); if(verb) raxa(az);
	iaxa(Fi,&ax,2); if(verb) raxa(ax);

	oaxa(Fo,&at,1);
	oaxa(Fo,&ag,2);
    }

    mapCC=sf_floatalloc2  (az.n,ax.n); LOOPCC( mapCC[ix][iz]=0.; );
    mapRC=sf_floatalloc2  (at.n,ag.n); LOOPRC( mapRC[ig][it]=0.; );
    rays =sf_complexalloc2(at.n,ag.n);
    
    sf_complexread(rays[0],at.n*ag.n,Fr);
    if(adj) {
	sf_floatread(mapRC[0],at.n*ag.n,Fi);
    } else {
	sf_floatread(mapCC[0],az.n*ax.n,Fi);
    }

    if(linear) {
	xmin = ax.o;
	xmax = ax.o + (ax.n-1)*ax.d;
	zmin = az.o;
	zmax = az.o + (az.n-1)*az.d;
	
	for(ig=0;ig<ag.n;ig++) {
	    if( ig%100 == 0 ) sf_warning("LINT %d of %d",ig,ag.n);
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
		    mapCC[ix][iz] += w00 * mapRC[ig][it];
		    mapCC[ix][jz] += w10 * mapRC[ig][it];
		    mapCC[jx][iz] += w01 * mapRC[ig][it];
		    mapCC[jx][jz] += w11 * mapRC[ig][it];
		} else {
		    mapRC[ig][it] += w00 * mapCC[ix][iz];
		    mapRC[ig][it] += w10 * mapCC[ix][jz];
		    mapRC[ig][it] += w01 * mapCC[jx][iz];
		    mapRC[ig][it] += w11 * mapCC[jx][jz];
		}
		
	    } /* it */
	} /* ig */
    } else {
	dr = sqrtf(az.d*az.d + ax.d*ax.d);

	if(! sf_getint("nsz",&nsz)) nsz=1; if(verb) sf_warning("nsz=%d",nsz);
	if(! sf_getint("nsx",&nsx)) nsx=1; if(verb) sf_warning("nsx=%d",nsx);

	/* loop over RC */
	for(ig=0;ig<ag.n;ig++) {
	    if( ig%100 == 0 ) sf_warning("SINT %d of %d",ig,ag.n);
	    for(it=0;it<at.n;it++) {
	
		z = cimagf(rays[ig][it]);
		x = crealf(rays[ig][it]);
		
		iz=floor((z-az.o)/az.d);
		ix=floor((x-ax.o)/ax.d);

		if(iz>=nsz && iz<=az.n-nsz-1 && ix>=nsx && ix<=ax.n-nsx-1) {
		    
		    /* sinc */
		    for(jz=iz-nsz; jz<=iz+nsz; jz++) {
			zo=az.o+jz*az.d;
			dz=zo-z;
			
			for(jx=ix-nsx; jx<=ix+nsx; jx++) {
			    xo=ax.o+jx*ax.d;
			    dx=xo-x;
			    
			    r=sqrt(dz*dz + dx*dx);
			    
			    if(SF_ABS(r) < 0.00001*dr) {
				sincr = 1.0;
			    } else {
				r = r / dr;
				sincr = sin(r)/r;
			    }

			    if(adj) {
				mapCC[jx][jz] += mapRC[ig][it] * sincr;
			    } else {
				mapRC[ig][it] += mapCC[jx][jz] * sincr;
			    }
			}
		    }
		    /* sinc end */

		} /* if iz,ix in bounds */

	    } /* it */
	} /* ig */
	
    } /* end linear */

    if(adj) {
	sf_floatwrite(mapCC[0],az.n*ax.n,Fo);
    } else {
	sf_floatwrite(mapRC[0],at.n*ag.n,Fo);
    }

}
