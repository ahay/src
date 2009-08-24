/* 1-D finite-difference wave extrapolation */
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
#include <math.h>
int main(int argc, char* argv[]) 
{
    int nx, ix;
    float dx;
    float *sig, *fsig, *v; 
    sf_file inp,vel, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    vel = sf_input("vel");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1= in input");

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);

    sig = sf_floatalloc(nx);
    fsig = sf_floatalloc(nx);
    v = sf_floatalloc(nx);
    sf_floatread(sig,nx,inp);		
    sf_floatread(v,nx,vel);		
    
     fsig[0] = sig[0]*2.0-sig[1]; 
     for (ix=1; ix < nx-1; ix++) {
          fsig[ix] = 2.0*sig[ix]-sig[ix-1]-sig[ix+1];
         }  
     
     fsig[nx-1] = sig[nx-1]*2.0-sig[nx-2]; 
     for (ix=0; ix < nx; ix++) {
          fsig[ix] = -fsig[ix]*v[ix]*v[ix];
         } 
    sf_floatwrite(fsig,nx,out);

    exit(0); 
}           
           
