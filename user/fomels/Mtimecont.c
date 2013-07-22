/* Forward or reverse time continuation using fast marching. */
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

#include<rsf.h>

#include "timecont.h"

int main (int argc,char* argv[]) 
{
    int n1, n2, n3, i, order,n123, n23, *p;
    float d1, d2, d3, slow;
    float *t, *t0, *v;
    bool isvel, forwd;
    sf_file vel, time, time0;

    sf_init (argc, argv);
    vel = sf_input("in");
    time = sf_output("out");
    time0 = sf_input("surf");

    if (SF_FLOAT != sf_gettype(vel)) 
	sf_error("Need float input");
    if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in input");
    if(!sf_histint(vel,"n3",&n3)) n3=1;

    if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in input");
    if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in input");
    if(!sf_histfloat(vel,"d3",&d3)) d3=d2;

    if(!sf_getbool("forwd",&forwd)) forwd=false;
    /* continue forward or backward */

    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* if y, the input is velocity; n, slowness squared */

    if(!sf_getint("order",&order)) order=2;
    /* [1,2] Accuracy order */

    n23 = n2*n3;
    n123 = n1*n23;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);
    t0 = sf_floatalloc (n23);
    p  = sf_intalloc   (n123);

    sf_floatread(v,n123,vel);
   /* transform velocity to slowness squared */
    if (isvel) {
	for(i = 0; i < n123; i++) {
	    slow = v[i];
	    v[i] = 1./(slow*slow);
	}
    } 
    
    sf_floatread(t0,n23,time0);

    timecont(t,t0,v,p,
	     n3,n2,n1,
	     d3,d2,d1,
	     order,forwd);
	
    sf_floatwrite (t,n123,time);

    exit (0);
}
