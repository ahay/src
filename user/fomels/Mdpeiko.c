/* 2-D eikonal solver based on dynamic programming. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "dynprog.h"

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, is;
    float d1, d2;
    float **vv, **tt;
    sf_file vel, ttime;

    sf_init(argc,argv);
    vel = sf_input("in");
    ttime = sf_output("out");

    if (!sf_histint(vel,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_getint("ishot",&is)) is=(n1-1)/2;
    /* shot location */

    vv = sf_floatalloc2(n1,n2);
    
    sf_floatread(vv[0],n1*n2,vel);    

    /* convert velocity to slowness */
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vv[i2][i1] = d1/vv[i2][i1];
	}
    }

    tt = dynprog_init(n2,n1,n1,d2/d1,true);

    dynprog(is,vv);

    sf_floatwrite(tt[0],n1*n2,ttime);

    exit(0);
}
    
