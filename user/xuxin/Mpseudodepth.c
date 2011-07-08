/* pseudo-depth mapping  */
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
#include <stdio.h>
#include "interp.h"

int main(int argc, char* argv[])
{
    sf_file Fi,Fo,Fh;  
    sf_axis ax,az,ah;
    int nx,nz,nh,ix,n;
    float oz,oh,dz,dh,**aa,**bb,**hh,*h;
    lint1dp L;

    sf_init(argc,argv);
    Fi = sf_input("in");    /* f(z) */
    Fh = sf_input("depth"); /* h(z) */
    Fo = sf_output("out");  /* f(h) */

    az = sf_iaxa(Fi,1); 
    ax = sf_iaxa(Fi,2); 
    nx = sf_n(ax); 
    nz = sf_n(az);
    dz = sf_d(az);
    oz = sf_o(az);
    if (!sf_histint(Fh,"n1",&n) || n != nz) sf_error("Need n1=%d in depth",nz);
    if (!sf_histint(Fh,"n2",&n) || n != nx) sf_error("Need n2=%d in depth",nx);

    if (!sf_getint("n",&nh))   nh = nz;
    /* target nz */
    if (!sf_getfloat("o",&oh)) oh = oz;
    /* target oz */
    if (!sf_getfloat("d",&dh)) dh = dz;
    /* target dz */    
    ah =sf_maxa(nh,oh,dh);

    aa = sf_floatalloc2(nz,nx); 
    bb = sf_floatalloc2(nh,nx);
    hh = sf_floatalloc2(nz,nx);

    sf_floatread(aa[0],nx*nz,Fi);  
    sf_floatread(hh[0],nx*nz,Fh);

    /* regular axis */
    h  = sf_floatalloc(nh);
    for (h[0]=oh, ix=1; ix < nh; ix++) {
	h[ix] = h[ix-1] + dh;
    }

    /* interpolation */
    for (ix=0; ix < nx; ix++) {
	L = lint1d_init(hh[ix],nz,h,nh);
	lint1d_extract(aa[ix],bb[ix],L);
	free(L);
    }

    /* write */
    sf_oaxa(Fo,ah,1);
    sf_oaxa(Fo,ax,2);
    sf_floatwrite(bb[0],nx*nh,Fo);

    exit (0);
}

