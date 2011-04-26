/* Diffraction focusing test. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "velcon.h"
#include "kurtosis.h"

int main(int argc, char* argv[]) 
{
    int nt, nx, ntx, n2, n3, next;
    float v0, v1, dt, dx, t0, kur;
    float **data;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) t0=0.;  

    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");

    ntx = nt*nx;

    if (!sf_getfloat("v0",&v0)) v0=SF_EPS;
    /* initial velocity */

    if (!sf_getfloat("v",&v1)) sf_error("Need v=");
    /* final velocity */

    if (!sf_getint("pad",&n2)) n2=nt; /* padding for stretch */
    if (!sf_getint("pad2",&n3)) n3=2*kiss_fft_next_fast_size((n2+1)/2);
    /* padding for FFT */

    if (!sf_getint("extend",&next)) next=4;
    /* trace extension */

    velcon_init(nt,nx,dt,dx,t0,n2,n3,next);

    data = sf_floatalloc2(nt,nx);

    sf_floatread(data[0],ntx,inp);

    kur = kurtosis(ntx,data[0]);
    sf_warning("kurtosis before: %g",kur);

    velcon(data,v0,v1);

    kur = kurtosis(ntx,data[0]);
    sf_warning("kurtosis after: %g",kur);

    sf_floatwrite(data[0],ntx,out);

    exit(0);
}
