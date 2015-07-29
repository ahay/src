/* Zero crossings. */
 
/*
  Copyright (C) 2011 King Abdullah University of Technology and Science, 
  Thuwal, Saudi Arabia.
 
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

int main(int argc, char* argv[])
{
    int it,nt, i,n2, levels;
    float *dat, *zc, d0,d1; /* dt, t0; */
	
    sf_file in, out; 
    
    /* initialization */
    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getint("levels",&levels)) levels = 3;  
    /*levels of quantization [2,3,5].
      levels=2	1: zero crossing or zero; 0: otherwise
      levels=3	1: positive to negative zc; -1 negative to positive zc; 0: otherwise
      levels=5	+/-2: positive/negative values; +/-1: as in levels=3; 0: zero. 
    */
    
    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1=");
    
    n2 = sf_leftsize(in,1);
	    
    dat = sf_floatalloc(nt);
    zc  = sf_floatalloc(nt);
		
    for (i=0; i<n2; i++) {
		
	sf_floatread(dat,nt,in);
	d1 = dat[0];
	zc[0] = 0.;
	for (it=1; it<nt; it++) {
	    d0 = d1;
	    d1 = dat[it];
	    if (levels==3) {
		if (d0*d1>=0.) {
		    zc[it] = 0.;
		} else {
		    if (d1>0.) zc[it] =  1.;
		    else 	   zc[it] = -1.;
		}
	    } else if (levels==2) {
		if (d0*d1>0.) zc[it] = 0.;
		else  zc[it] =  1.;
	    } else {
		if (d0*d1>0.) {
		    if (d1>0.) zc[it] =  2.;
		    else 	   zc[it] = -2.;
		} else if (d0*d1<0.) { 
		    if (d1>0.) zc[it] =  1.;
		    else	   zc[it] = -1.;
		} else {
		    zc[it] = 0.;
		}
	    }
	} /* nt */
		
	sf_floatwrite(zc,nt,out);
    }
    exit(0);
}



