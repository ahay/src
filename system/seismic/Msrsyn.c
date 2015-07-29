/* Synthesize shot/receiver wavefields for 3-D SR migration */
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

#include <rsf.h>

int main (int argc, char *argv[])
{
    sf_axis aw;
    sf_axis arx,asx,ax;
    sf_axis ary,asy,ay;
    sf_axis ae;

    sf_file Fs=NULL;  /* source   file sou[nw]         */
    sf_file Fr=NULL;  /* receiver file rec[nw][nrx][nry][nsx][nsy] */
    sf_file Fsw=NULL; /* source wavefield [nw][ nx][ ny][ne] */
    sf_file Frw=NULL; /* source wavefield [nw][ nx][ ny][ne] */

    sf_complex   *rr=NULL;
    sf_complex   *ss=NULL;
    sf_complex ***rw=NULL;
    sf_complex ***sw=NULL;
    int   isx,isy,irx,iry,ix,iy,iw;
    int   nsx,nsy,nrx,nry,nx,ny,nw;
    float ox,dx,oy,dy;
    float  sx, sy, rx, ry,xx,yy;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    Fr = sf_input ( "in");
    Fs = sf_input ("wav");
    Fsw= sf_output("swf"); sf_settype(Fsw,SF_COMPLEX);
    Frw= sf_output("out"); sf_settype(Frw,SF_COMPLEX);

    if (!sf_getint  ("nx",&nx)) sf_error ("Need nx="); /* x samples */
    if (!sf_getfloat("dx",&dx)) sf_error ("Need dx="); /* x sampling */
    if (!sf_getfloat("ox",&ox)) sf_error ("Need ox="); /* x origin */
    ax = sf_maxa(nx,ox,dx); sf_setlabel(ax,"x");

    if (!sf_getint  ("ny",&ny)) ny=1; /* y samples */
    if (!sf_getfloat("dy",&dy)) dy=1; /* y sampling */
    if (!sf_getfloat("oy",&oy)) oy=0; /* y origin */
    ay = sf_maxa(ny,oy,dy); sf_setlabel(ay,"y");

    aw  = sf_iaxa(Fr,1); nw =sf_n(aw) ; sf_setlabel(aw, "w");
    arx = sf_iaxa(Fr,2); nrx=sf_n(arx); sf_setlabel(arx,"rx");
    ary = sf_iaxa(Fr,3); nry=sf_n(ary); sf_setlabel(ary,"ry");
    asx = sf_iaxa(Fr,4); nsx=sf_n(asx); sf_setlabel(asx,"sx");
    asy = sf_iaxa(Fr,5); nsy=sf_n(asy); sf_setlabel(asy,"sy");

    /* experiments axis */
    ae = sf_maxa(nsx*nsy,0,1); sf_setlabel(ae,"e");

    sf_oaxa(Fsw, aw,1); sf_oaxa(Frw, aw,1);
    sf_oaxa(Fsw, ax,2); sf_oaxa(Frw, ax,2);
    sf_oaxa(Fsw, ay,3); sf_oaxa(Frw, ay,3);
    sf_oaxa(Fsw, ae,4); sf_oaxa(Frw, ae,4);

    /*------------------------------------------------------------*/
    ss = sf_complexalloc (nw);
    rr = sf_complexalloc (nw);
    sw = sf_complexalloc3(nw,nx,ny);
    rw = sf_complexalloc3(nw,nx,ny);

    /* SOURCE wavefield */
    sf_complexread(ss,nw,Fs);

    for(    isy=0;isy<nsy;isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	for(isx=0;isx<nsx;isx++) { sx = sf_o(asx) + isx * sf_d(asx);

	    for        (iy=0;iy<ny;iy++) {
		for    (ix=0;ix<nx;ix++) {
		    for(iw=0;iw<nw;iw++) {
			sw[iy][ix][iw] = sf_cmplx(0.,0.);
		    }
		}
	    }

	    yy = sy; iy = 0.5+(yy-oy)/dy;
	    xx = sx; ix = 0.5+(xx-ox)/dx;
	    if(iy>=0 && iy<ny && ix>=0 && ix<nx) {
		
		for( iw=0;iw<nw;iw++) {
		    sw[iy][ix][iw] = ss[iw];
		}
	    }

	    sf_complexwrite(sw[0][0],nx*ny*nw,Fsw);
	}
    }

    /* RECEIVER wavefield */
    for(    isy=0;isy<nsy;isy++) { sy = sf_o(asy) + isy * sf_d(asy);
	for(isx=0;isx<nsx;isx++) { sx = sf_o(asx) + isx * sf_d(asx);

	    for        (iy=0;iy<ny;iy++) {
		for    (ix=0;ix<nx;ix++) {
		    for(iw=0;iw<nw;iw++) {
			rw[iy][ix][iw] = sf_cmplx(0.,0.);
		    }
		}
	    }

	    for(    iry=0;iry<nry;iry++) { ry = sf_o(ary) + iry * sf_d(ary);
		for(irx=0;irx<nrx;irx++) { rx = sf_o(arx) + irx * sf_d(arx);

		    sf_complexread(rr,nw,Fr);

		    yy = sy + ry; iy = 0.5+(yy-oy)/dy;
		    xx = sx + rx; ix = 0.5+(xx-ox)/dx;

		    if(iy>=0 && iy<ny && ix>=0 && ix<nx) {
			for( iw=0;iw<nw;iw++) {
			    rw[iy][ix][iw] = rr[iw];
			}
		    }
		}
	    }
	    sf_complexwrite(rw[0][0],nx*ny*nw,Frw);
	}
    }
    /*------------------------------------------------------------*/

    ;         ;            free(ss);
    ;         ;            free(rr);
    free(**sw); free(*sw); free(sw);
    free(**rw); free(*rw); free(rw);

    exit (0);
}
