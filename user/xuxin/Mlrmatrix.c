/* Compute fullrank/lowrank extrapolator matrix */
/* This only works for SMALL size models, otherwise it fails to allocate memory. */
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
#include <math.h>

int main(int argc, char* argv[])
{
    int i,j,m,n,nx,nz,nkx,nkz;
    float dt,okx,okz,dkx,dkz,kx,kz,*v,*k,*w;
    sf_file Fi,Fk,Fo,Fl,Fm,Fr;
    sf_axis ax,az,akx,akz;

    sf_init(argc,argv);

    if (!sf_getfloat("dt",&dt)) dt=1.0;
    /* time step */

    Fi = sf_input("in"); /* velocity */
    Fk = sf_input("fft");
    Fo = sf_output("out");
/*
    Fl = sf_input("left");
    Fm = sf_input("middle");
    Fr = sf_input("right");
*/
    
    /* read wavenubmers */
    akx= sf_iaxa(Fk,1); 
    akz= sf_iaxa(Fk,2); 
    nkx= sf_n(akx);   nkz= sf_n(akz);
    okx= sf_o(akx);   okz= sf_o(akz);
    dkx= sf_d(akx);   dkz= sf_d(akz);
    m = nkx*nkz;
/*
    kx= sf_floatalloc(nkx);
    kz= sf_floatalloc(nkz);
    for (ix=0; ix < nkx; ix++) {
	kx[ix] = okx + ix*dkx;
    }
    for (iz=0; iz < nkz; iz++) {
	kz[iz] = okz + iz*dkz;
    }
*/
    k = sf_floatalloc(m);
    for (i=0; i < nkz; i++) {
	kz = okz + i*dkz;
	for (j=0; j < nkx; j++) {
	    kx = okx + j*dkx;
	    k[j+i*nkx] = 2*SF_PI*hypot(kx,kz);
	}
    }

    /* read velocity */
    ax = sf_iaxa(Fi,1);    nx = sf_n(ax);   
    az = sf_iaxa(Fi,2);    nz = sf_n(az);
    n = nx*nz;
    v = sf_floatalloc(n);
    sf_floatread(v,n,Fi);

    /* fullrank matrix */
    w = sf_floatalloc(m*n);
    for (i=0; i < n; i++) {
	for (j=0; j < m; j++) {
	    w[j+i*m] = 2.0f*cos(v[i]*k[j]*dt);
	}
    }


    /* n1 = nx*nz; */
    /* n2 = nkx*nkz; */
    Fl = sf_input("left");
    Fr = sf_input("right");
    Fm = sf_input("middle");

    if (!sf_histint(Fm,"n1",&m1)) sf_error("No n1= in middle");
    if (!sf_histint(Fm,"n2",&m2)) sf_error("No n2= in middle");
    if (!sf_histint(Fl,"n1",&n) || n != nx*nz) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(Fl,"n2",&n) || n != m1) sf_error("Need n2=%d in left",m1);
    if (!sf_histint(Fr,"n1",&n) || n != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Fr,"n2",&n) || n != n2) sf_error("Need n2=%d in right",nk);
 
    ll = sf_floatalloc2(n1,m1);
    mm = sf_floatalloc2(m1,m2);
    rr = sf_floatalloc2(m2,n2);
    
    sf_floatread(ll[0],n1*m1,Fl);
    sf_floatread(mm[0],m1*m2,Fm);
    sf_floatread(rr[0],m2*n2,Fr);
    
    /* mm*rr */
    mmrr = sf_floatalloc2(m1,n2);
    for (im=0; im < m1; im++) {
	for (ik=0; ik < n2; ik++) {
	    for (s=0, i=0; i < m2; i++)
		s += mm[i][im]*rr[ik][i];
	    mmrr[ik][im] = s;
	}
    }
    /* lowrank matrix = ll*mm*rr */
    aa = sf_floatalloc2(n1,n2);
    for (i=0; i < n1; i++) {
	for (j=0; j < n2; j++) {
	    for (s=0, im=0; im < m1; im++)
		s += ll[im][i]*mmrr[j][im];
	    aa[j][i] = s;
	}
    }
    
    /* write */
    sf_putint(Fo,"n1",m);
    sf_putint(Fo,"n2",n);
    sf_putstring(Fo,"label1","wavenumber");
    sf_putstring(Fo,"label2","space");
    sf_putstring(Fo,"unit1","1/");
    sf_putstring(Fo,"unit2","");
    sf_floatwrite(w,m*n,Fo);
}
 
