/* Synthesize shot/receiver wavefields for 3-D SR migration */
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

#define INDEX(x,a) 0.5+(x-a.o)/a.d;

int main (int argc, char *argv[])
{
    axa aw;
    axa arx,asx,ax;
    axa ary,asy,ay;
    axa ae;

    sf_file Fs;  /* source   file sou[nw]         */
    sf_file Fr;  /* receiver file rec[nw][nrx][nry][nsx][nsy] */
    sf_file Fsw; /* source wavefield [nw][ nx][ ny][ne] */
    sf_file Frw; /* source wavefield [nw][ nx][ ny][ne] */

    float complex   *rr;
    float complex   *ss;
    float complex ***rw;
    float complex ***sw;
    int   isx,isy,irx,iry,ix,iy,iw;
    float  sx, sy, rx, ry,xx,yy;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    Fr = sf_input ( "in");
    Fs = sf_input ("wav");
    Fsw= sf_output("swf"); sf_settype(Fsw,SF_COMPLEX);
    Frw= sf_output("out"); sf_settype(Frw,SF_COMPLEX);
    
    if (!sf_getint  ("nx",&ax.n)) sf_error ("Need nx=");
    if (!sf_getfloat("dx",&ax.d)) sf_error ("Need dx=");
    if (!sf_getfloat("ox",&ax.o)) sf_error ("Need ox=");
    ax.l="x";

    if (!sf_getint  ("ny",&ay.n)) ay.n=1;
    if (!sf_getfloat("dy",&ay.d)) ay.o=0;
    if (!sf_getfloat("oy",&ay.o)) ay.d=1;
    ay.l="y";

    iaxa(Fr,&aw ,1);  aw.l= "w";
    iaxa(Fr,&arx,2); arx.l="rx";
    iaxa(Fr,&ary,3); ary.l="ry";
    iaxa(Fr,&asx,4); asx.l="sx";
    iaxa(Fr,&asy,5); asy.l="sy";

    /* experiments axis */
    ae.o=0;
    ae.d=1;
    ae.n=asx.n*asy.n;
    ae.l="e";

    oaxa(Fsw, &aw,1); oaxa(Frw, &aw,1);
    oaxa(Fsw, &ax,2); oaxa(Frw, &ax,2);
    oaxa(Fsw, &ay,3); oaxa(Frw, &ay,3);
    oaxa(Fsw, &ae,4); oaxa(Frw, &ae,4);

    sf_close();

/*------------------------------------------------------------*/

    ss = sf_complexalloc (aw.n);
    rr = sf_complexalloc (aw.n);
    sw = sf_complexalloc3(aw.n,ax.n,ay.n);
    rw = sf_complexalloc3(aw.n,ax.n,ay.n);

    /* SOURCE wavefield */
    sf_complexread(ss,aw.n,Fs);

    for( isy=0;isy<asy.n;isy++) {             
	sy = asy.o + isy * asy.d;
	for( isx=0;isx<asx.n;isx++) {         
	    sx = asx.o + isx * asx.d;
	    
	    for(iy=0;iy<ay.n;iy++) {
		for(ix=0;ix<ax.n;ix++) {
		    for(iw=0;iw<aw.n;iw++) {
			sw[iy][ix][iw] = 0.;
		    }
		}
	    }
	    
	    yy = sy; iy = INDEX(yy,ay);
	    xx = sx; ix = INDEX(xx,ax);
	    if(iy>=0 && iy<ay.n && ix>=0 && ix<ax.n) {
		
		for( iw=0;iw<aw.n;iw++) {
		    sw[iy][ix][iw] = ss[iw];
		}
	    }
	    
	    sf_complexwrite(sw[0][0],ax.n*ay.n*aw.n,Fsw);
	}
    }

    /* RECEIVER wavefield */
    for( isy=0;isy<asy.n;isy++) {	      
	sy = asy.o + isy * asy.d;
	for( isx=0;isx<asx.n;isx++) {         
	    sx = asx.o + isx * asx.d;
	    
	    for(iy=0;iy<ay.n;iy++) {
		for(ix=0;ix<ax.n;ix++) {
		    for(iw=0;iw<aw.n;iw++) {
			rw[iy][ix][iw] = 0.;
		    }
		}
	    }
	    
	    for( iry=0;iry<ary.n;iry++) {     
		ry = ary.o + iry * ary.d;
		for( irx=0;irx<arx.n;irx++) { 
		    rx = arx.o + irx * arx.d;

		    sf_complexread(rr,aw.n,Fr);
		    
		    yy = sy + ry; iy = INDEX(yy,ay);
		    xx = sx + rx; ix = INDEX(xx,ax);
		    
		    if(iy>=0 && iy<ay.n && ix>=0 && ix<ax.n) {
			for( iw=0;iw<aw.n;iw++) {
			    rw[iy][ix][iw] = rr[iw];
			}
		    }
		}
	    }
	    sf_complexwrite(rw[0][0],ax.n*ay.n*aw.n,Frw);
	}
    }
/*------------------------------------------------------------*/
    
    ;         ;            free(ss);
    ;         ;            free(rr);
    free(**sw); free(*sw); free(sw);
    free(**rw); free(*rw); free(rw);
    
    exit (0);
}
