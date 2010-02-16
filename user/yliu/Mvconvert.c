/* 2-D velocity mapping from manual picking to rsf RMS format. 
Covert from
--------------------------------------
Time_i(ms) RMS_i(m/s) Lateral_j(trace)
...        ...        ...
-1         0          Lateral_j
...        ...        ...
--------------------------------------
       to
regular RSF_RMS velocity grid

Specify either n1= o1= d1= n2= o2= d2= or pattern=
*/
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
    
    int pn1, pn2, pn3, n1, n2, i1, i2, pi2, pi3;
    int index1, index2, begt, endt, begx, endx, bt, bx, shift;
    
    float *trace, *vel, d1, d2, o1, o2;
    sf_file in, out, pattern=NULL;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&pn1)) sf_error("No pn1= in input");
    if (!sf_histint(in,"n2",&pn2)) sf_error("No pn2= in input");
    pn3 = sf_leftsize(in,2);

    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }

    if (!sf_getint("n1",&n1) && 
	(NULL== pattern ||
	 !sf_histint(pattern,"n1",&n1))) sf_error("Need n1=");
    /* Output grid size n1 */
    if (!sf_getfloat("d1",&d1) && 
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d1",&d1))) sf_error("Need d1=");
    /* Output sampling d1 */
    if (d1 > 0.1) d1/=1000.;
    if (!sf_getfloat("o1",&o1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o1",&o1))) sf_error("Need o1=");
    /* Output origin o1 */
    if (!sf_getint("n2",&n2) && 
	(NULL== pattern ||
	 !sf_histint(pattern,"n2",&n2))) sf_error("Need n2=");
    /* Output grid size n2 */
    if (!sf_getfloat("d2",&d2) && 
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d2",&d2))) sf_error("Need d2=");
    /* Output sampling d2 */
    if (!sf_getfloat("o2",&o2) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o2",&o2))) sf_error("Need o2=");
    /* Output origin o2 */

    if (!sf_getint("shift",&shift)) shift=0;
    /* Lateral shift */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"o2",o2);
    sf_putfloat(out,"d2",d2);
    sf_putstring(out,"label1","Time");
    sf_putstring(out,"unit1","");
    sf_putstring(out,"label2","Midpoint");
    sf_putstring(out,"unit2","");
       
    trace = sf_floatalloc(pn1);
    vel = sf_floatalloc(n1*n2);
    index1 = 0;
    index2 = 0;
    endx = -1;
    for(pi3=0; pi3 < pn3; pi3++) {
	begt = 0;
	begx = 0;
        for(pi2=0; pi2 < pn2; pi2++) {
	    sf_floatread(trace,pn1,in);
	    if(-1==trace[0]) {
		begx = endx+1;
		begt = 0;
	    } else {
		bt = (int) (trace[0]/(d1*1000.)+0.5);
		endt = (n1 <= bt)? n1:bt;
		bx = (int) (trace[2]+0.5);
		endx = (n2 <= bx)? n2:(bx+shift);
		if (index1 <= endt) index1 = endt;
		if (index2 <= endx) index2 = endx;
		for (i2=begx; i2 <= endx; i2++) {
		    for (i1=begt; i1 <= endt; i1++) {
			vel[i2*n1+i1] = trace[1]/1000.;
		    }
		}
		begt = endt+1;
	    }
	}
	if ((n2-1) > index2) {
	    for (i2=index2; i2 < n2; i2++) {
	        for (i1=0; i1 < n1; i1++) {
		  vel[i2*n1+i1] = vel[(index2-1)*n1+i1];
		}
	    }
	}
	if ((n1-1) > index1) {
	    for (i2=0; i2 < n2; i2++) {
	        for (i1=index1; i1 < n1; i1++) {
		    vel[i2*n1+i1] = vel[i2*n1+index1-1];
		}
	    }
	}
	    
	sf_floatwrite(vel,n1*n2,out);
    }

    exit (0);
}

/* 	$Id$	 */
