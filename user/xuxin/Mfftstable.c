/* stablization trick */
/* This code has problem. Sometimes when use computed left as input again will give beta still greater than 2, by 0.0001. Either the computation is wrong or double precision is required. */
/*
  Copyright (C) 2011 KAUST
  
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
#include "fft2.h"

int main(int argc, char* argv[])
{
    sf_file Fo,Fl,Fm,Fr,Fb;
    float **ll,**mm,**rr,**lm,*b,c,M;
    int i,j,n,m1,m2,nx,nz,nx2,nz2,nkx,nkz,nn,nk,pad1,gd;
    bool cmplx;

    sf_init(argc,argv);

    Fm = sf_input("in");
    Fl = sf_input("left");
    Fr = sf_input("right");
    Fo = sf_output("out");
    Fb = sf_output("b");

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false;
    /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */
    if (!sf_getint("nx",&nx)) sf_error("No nx= in input");
    /* nx */
    if (!sf_getint("nz",&nz)) sf_error("No nz= in input");
    /* nz */
    if (!sf_getfloat("M",&M)) M=2.0f;
    /* max abs */

    /* num of x */
    nn = nx*nz;

    /* num of k */
    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nkx = nx2;
    nkz = (int)(nk/nkx);
    if (nkx & 1 || nkz & 1) sf_error("Need nkx=even and nkz=even");

    /* lowrank matrices */
    if (!sf_histint(Fm,"n1",&m1)) sf_error("No n1= in middle");
    if (!sf_histint(Fm,"n2",&m2)) sf_error("No n2= in middle");
    if (!sf_histint(Fl,"n1",&n) || n != nn) sf_error("Need n1=%d in left",nn);
    if (!sf_histint(Fl,"n2",&n) || n != m1) sf_error("Need n2=%d in left",m1);
    if (!sf_histint(Fr,"n1",&n) || n != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Fr,"n2",&n) || n != nk) sf_error("Need n2=%d in right",nk); 

    ll = sf_floatalloc2(nn,m1);
    mm = sf_floatalloc2(m1,m2);
    rr = sf_floatalloc2(m2,nk);

    sf_floatread(ll[0],nn*m1,Fl);
    sf_floatread(mm[0],m1*m2,Fm);
    sf_floatread(rr[0],m2*nk,Fr);

    /* ll*mm */
    lm = sf_floatalloc2(nn,m2);
    for (i=0; i < nn; i++) {
	for (j=0; j< m2; j++) {
	    c = 0.0f;
	    for (n=0; n < m1; n++)
		c += ll[n][i]*mm[j][n];
	    lm[j][i] = c;
	}
    }

    /* lowrank at k=0 */
    b = sf_floatalloc(nn);
    gd = (int)(nkx/2)*nkz + (int)(nkz/2); /* k=0 index */
    for (i=0; i < nn; i++) {
	c = 0;
	for (j=0; j < m2; j++) {
	    c += lm[j][i]*rr[gd][j];
	}
	b[i] = c;
    }

    /* apply */
    for (i=0; i < nn; i++) {
	for (j=0; j < m1; j++) {
	    c = abs(b[i]);
	    ll[j][i] *= (c > M) ? M/c : 1.0;
	}
    }

    /* write updated left */
    sf_putint(Fo,"n1",nn);
    sf_putint(Fo,"n2",m1);
    sf_floatwrite(ll[0],nn*m1,Fo);

    /* write b */
    sf_putint(Fb,"n1",nz);
    sf_putint(Fb,"n2",nx);
    sf_putstring(Fb,"label1","z");
    sf_putstring(Fb,"label2","x");
    sf_floatwrite(b,nn,Fb);

    exit(0);
}
