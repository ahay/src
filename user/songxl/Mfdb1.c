/* 1-D Finite-difference wave extrapolation */
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
#include <rsf.h>
#include "abcpass.h"
#include "lap.h"

int main(int argc, char* argv[]) 
{
    int nx, nt, ix, it, nb, nxb, abc, order;
    float dt, dx;
    float *old, *nxt, *cur, *sig, *v, 
    *w, *dercur, *derold;
    sf_file in, out, vel;

    sf_init(argc,argv);
    in  = sf_input("in");
    vel = sf_input("vel");   /* velocity */
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("nb",&nb)) nb =20; 
    if (!sf_getint("abc",&abc)) abc =0; /*absorbing boundary condition 1: cos 0: exp*/ 
    if (!sf_getint("order",&order)) order=2; /*FD order: 2,4*/ 

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
//    sf_putfloat(out,"o1",x0);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 

    nxb=nx+2*nb;


    sig = sf_floatalloc(nx);
    old = sf_floatalloc(nxb);
    nxt = sf_floatalloc(nxb);
    cur = sf_floatalloc(nxb);
    dercur = sf_floatalloc(nxb);
    derold = sf_floatalloc(nxb);
    v = sf_floatalloc(nxb);

    w = sf_floatalloc(nb);
   /* sb = 4.0*nb;
    for(ib=0; ib<nb; ib++){
       fb = ib/(sqrt(2.0)*sb);
       w[ib] = exp(-fb*fb);
    }*/
    abc_cal(abc,nb,0.0001,w);

    sf_floatread(v,nx,vel);
    sf_floatread(sig,nx,in);		
    sf_floatwrite(sig,nx,out);

    for (ix= nb+nx-1; ix > nb-1; ix--) {
        v[ix] = v[ix-nb];
         } 
    
    for (ix=0; ix < nb; ix++){
        v[ix] = v[nb];
        v[ix+nb+nx] = v[nx+nb-1];
        }


       /* for (ix=0; ix < nb; ix++){
            aa[ix][0] *= (1.0+cosf(((float)(nb-1.0)-(float)ix)/(float)(nb-1)*pi))/2.0;
            aa[ix+nb+nx][1] *= (1.0+cosf((float)ix/(float)(nb-1)*pi))/2.0;
        }*/
	/* initial conditions */
    for (ix=0; ix < nx; ix++){
        cur[ix+nb] =  sig[ix];
    }
    for (ix=0; ix < nb; ix++){
        cur[ix] = cur[nb];
        cur[ix+nb+nx] = cur[nx+nb-1];
        }
    for (ix=0; ix < nxb; ix++){
        old[ix] =  0.0; 
	nxt[ix] = 0.;
        derold[ix] = cur[ix]/dt ;
    }

    free(sig);
    lap1_init(nxb);

    /* propagation in time */
    for (it=1; it < nt; it++) {
/*	nxt[0] = cur[1]-2*cur[0];
	for (ix=1; ix < nxb-1; ix++) {
	    nxt[ix] = cur[ix+1]+cur[ix-1]-2.0*cur[ix]; 
	}
	nxt[nxb-1] = cur[nxb-2]-2.0*cur[nxb-1]; */
        lap1(cur,nxt,order);
	for (ix=0; ix < nxb; ix++){
            dercur[ix]= derold[ix] + nxt[ix]*v[ix]*v[ix]*dt/(dx*dx);
            } 
	for (ix=0; ix < nxb; ix++) {
	    nxt[ix] =  cur[ix] + dercur[ix]*dt;
	}
        for (ix=0; ix < nb; ix++){
            nxt[ix] *= w[ix];
            nxt[ix+nb+nx] *= w[nb-1-ix];
            dercur[ix] *= w[ix];
            dercur[ix+nb+nx] *= w[nb-1-ix];
        }
	sf_floatwrite(nxt+nb,nx,out);
	for (ix=0; ix < nxb; ix++) {
	    old[ix] = cur[ix];
	    cur[ix] = nxt[ix];
	    derold[ix] = dercur[ix];
	}
    }

    sf_fileclose(vel);
    sf_fileclose(in);
    sf_fileclose(out);
    exit(0);
}
