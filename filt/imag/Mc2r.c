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

#define LOOPCC(a) for(iz=0;iz<az.n;iz++){ for(ix=0;ix<ax.n;ix++){  {a} }}
#define LOOPRC(a) for(it=0;it<at.n;it++){ for(ig=0;ig<ag.n;ig++){  {a} }}

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

    mapCC=sf_floatalloc2  (ax.n,az.n); LOOPCC( mapCC[iz][ix]=0.; );
    mapRC=sf_floatalloc2  (ag.n,at.n); LOOPRC( mapRC[it][ig]=0.; );

    rays =sf_complexalloc2(ag.n,at.n);
    sf_complexread(rays[0],ag.n*at.n,Fr);
    if(adj) {
	sf_floatread(mapRC[0],ag.n*at.n,Fi);
    } else {
	sf_floatread(mapCC[0],ax.n*az.n,Fi);
    }

    if(linear) {
	xmin = ax.o;
	xmax = ax.o + (ax.n-1)*ax.d;
	zmin = az.o;
	zmax = az.o + (az.n-1)*az.d;
	
	for(it=0;it<at.n;it++) {		
	    if( it%100 == 0 ) sf_warning("LINT %d of %d",it,at.n);	    

	    for(ig=0;ig<ag.n;ig++) {

		z = cimagf(rays[it][ig]);
		x = crealf(rays[it][ig]);
		
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
		  00 . 01    --> x
	          .    .    |
		  10 . 11   v
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
		    mapCC[iz][ix] += w00 * mapRC[it][ig];
		    mapCC[iz][jx] += w10 * mapRC[it][ig];
		    mapCC[jz][ix] += w01 * mapRC[it][ig];
		    mapCC[jz][jx] += w11 * mapRC[it][ig];
		} else {
		    mapRC[it][ig] += w00 * mapCC[iz][ix];
		    mapRC[it][ig] += w10 * mapCC[iz][jx];
		    mapRC[it][ig] += w01 * mapCC[jz][ix];
		    mapRC[it][ig] += w11 * mapCC[jz][jx];
		}
		
	    } /* ig */
	} /* it */
    } else {
	dr = sqrtf(az.d*az.d + ax.d*ax.d);

	if(! sf_getint("nsz",&nsz)) nsz=1; if(verb) sf_warning("nsz=%d",nsz);
	if(! sf_getint("nsx",&nsx)) nsx=1; if(verb) sf_warning("nsx=%d",nsx);

	/* loop over RC */
	for(it=0;it<at.n;it++) {
	    if( it%100 == 0 ) sf_warning("SINT %d of %d",it,at.n);
	    for(ig=0;ig<ag.n;ig++) {
	
		z = cimagf(rays[it][ig]);
		x = crealf(rays[it][ig]);
		
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
				mapCC[jz][jx] += mapRC[it][ig] * sincr;
			    } else {
				mapRC[it][ig] += mapCC[jz][jx] * sincr;
			    }
			}
		    }
		    /* sinc end */

		} /* if iz,ix in bounds */

	    } /* ig */
	} /* it */
	
    } /* end linear */

    if(adj) {
	sf_floatwrite(mapCC[0],ax.n*az.n,Fo);
    } else {
	sf_floatwrite(mapRC[0],ag.n*at.n,Fo);
    }

}
