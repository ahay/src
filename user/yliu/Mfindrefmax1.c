/* Find the sampled position of max value after reference point along fast dimension. */
/*
  Copyright (C) 2022 Jilin University
  
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

float find_max (const float* a, int n, int* max_ind_pt);

int main (int argc, char* argv[])
{
    bool verb;
    
    float o1, d1;	
    int   n1, n2, i2, i1, zp;
    int   maxpo=0;
    float max_value;
    float *column,*max; 
    int *zepo, *max_pt;	
    sf_file inp, out, max_val, zeropo;
	
    sf_init (argc,argv);
	
    inp = sf_input("in");
    out = sf_output("out");
    
    sf_settype (out, SF_INT);

    if (!sf_getbool("verb",&verb)) verb=false;
	
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");
    
    n2=sf_leftsize(inp,1);
    sf_unshiftdim(inp, out, 1);

    if (NULL != sf_getstring("max_val")) {
	max_val=sf_output("max_val");
    } else {
	max_val= NULL;
    }

    if (NULL != sf_getstring("zeropo")) {
	zeropo = sf_input("zeropo");
    } else {
	zeropo = NULL;
    }
    
    sf_putint(out,"n1",n2);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"d1",1);
    
    sf_putint(out,"n2",1);
    sf_putfloat(out,"o2",0.);
    sf_putfloat(out,"d2",1);

    column = sf_floatalloc(n1);
    max_pt = sf_intalloc(n2);
    max = sf_floatalloc(n2);
    zepo = sf_intalloc(n2);
	
    if (NULL!=max_val) {
	sf_putint(max_val,"n1",n2);
    	sf_putfloat(max_val,"o1",0.);
	sf_unshiftdim(inp, max_val, 1);
    }

    if(NULL != zeropo) {
	sf_intread(zepo,n2,zeropo);
    }

    for (i2=0; i2<n2; i2++) { 	
	if(verb) sf_warning("trace %d of %d;",i2+1,n2);
	
	sf_floatread(column,n1,inp);

	if(NULL != zeropo) {
	    zp = zepo[i2];
	    max_value = column[zp+1];
	    for (i1=zp+1; i1<n1; i1++) {
		if (column[i1] >= max_value) {
		    max_value = column[i1];
		    max[i2] = max_value;
		    maxpo = i1;
		}
	    }
	} else {
	    max_value = column[0];
	    for (i1=0; i1<n1; i1++) {
		if (column[i1] >= max_value) {
		    max_value = column[i1];
		    max[i2] = max_value;
		    maxpo = i1;
		}
	    }
	}
	max_pt[i2]  = maxpo;
    } 

    sf_intwrite(max_pt,n2,out);
    
    if (NULL!=max_val) {
	sf_floatwrite(max,n2,max_val);
    }
    
    exit (0);
}

