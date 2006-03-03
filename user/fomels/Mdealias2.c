/* 2-D trace interpolation to a denser grid using PWD.

It may be necessary to bandpass the data before and after dealiasing 
to ensure that the temporal spectrum is banded. Rule of thumb: if 
max(jx,jy)=N, the temporal bandwidth should be 1/N of Nyquist.
*/
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
#include "interpd.h"

int main (int argc, char *argv[])
{
    int n1, n2, n12, m2, m12;
    float d2, **u1, **uu1, **p;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if(SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    m2 = n2*2-1; 
    if (n2 > 1) d2 /= 2;
    
    n12 = n1*n2;
    m12 = n1*m2;
    
    sf_putint(out,"n2",m2); 
    sf_putfloat(out,"d2",d2);

    u1 = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    uu1 = sf_floatalloc2(n1,m2);

    interp_init (n1, 0.0001);

    sf_floatread(u1[0],n12,in);
    sf_floatread(p[0],n12,dip);

    interp2(n2,u1,uu1,p);
    
    sf_floatwrite(uu1[0],m12,out);
  
    exit (0);
}

/* 	$Id$	 */
