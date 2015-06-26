/* Simple cube fault synthetic */
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

static void remove_mean(int n, float *a)
{
    int i;
    float t;

    t = 0.;
    for (i=0; i < n; i++) {
	t += a[i];
    }
    t /= n;
    for (i=0; i < n; i++) {
	a[i] -= t;
    }
}

int main (int argc, char* argv[])
{
    int n, x, y, z, z1;
    float p, q, w;
    float *trace, *t1, *t2;
    sf_triangle tr;
    sf_file cube;

    sf_init (argc,argv);
    cube = sf_output("out");

    if (!sf_getint("n", &n)) n=51;
    /* cube dimensions */

    sf_setformat(cube,"native_float");
    sf_putint(cube,"n1",n);
    sf_putint(cube,"n2",n);
    sf_putint(cube,"n3",n);
    sf_putfloat(cube,"d1",1.0);
    sf_putfloat(cube,"d2",1.0);
    sf_putfloat(cube,"d3",1.0);
    sf_putfloat(cube,"o1",0.0);
    sf_putfloat(cube,"o2",0.0);
    sf_putfloat(cube,"o3",0.0);

    if (!sf_getfloat("p",&p)) p=0.5;
    /* inline slope */
    if (!sf_getfloat("q",&q)) q=0.5; 
    /* crossline slope */

    init_genrand(2008);

    trace = sf_floatalloc(n);
    t1 = sf_floatalloc(n);
    t2 = sf_floatalloc(3*n);

    sf_random (n,t1);
    sf_random (3*n,t2);
  
    tr = sf_triangle_init(2,n,false);
    sf_smooth2 (tr,0,1,false,t1);
    sf_triangle_close(tr);
    remove_mean(n,t1);

    tr = sf_triangle_init(2,3*n,false);
    sf_smooth2 (tr,0,1,false,t2);
    sf_triangle_close(tr);
    remove_mean(3*n,t2);
    
    for (y=0; y < n; y++) {
	for (x=0; x < n; x++) {
	    for (z=0; z < n; z++) {
		if (z < 3*n/2 - x - y - 3) {
		    trace[z] = t1[z];
		} else {
		    w = 2*n + z - p*(x+1) - q*(y+1); 
		    z1 = floorf(w);  w -= z1;
		    trace[z] = (1.-w) * t2[z1] + w * t2[z1+1];
		}
	    }
	    sf_floatwrite (trace,n,cube);
	}
    }

    exit(0);
}

