/* Recover all equivalent Q values according to reference point and non-zero point. */
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

#include <rsf.h>
#include <math.h>

int main(int argc, char *argv[])
{
    bool  verb;

    int   n1, n2, i1, i2;                    
    float d1, o1;
    int   spo, mspo, zp, *pos, *zepo;
    float *oeq, *aeq;
    sf_file inp, out, zpo, zeropo;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;

    if (!sf_histint(inp, "n1", &n1)) sf_error("No n1=in input");
    if (!sf_histfloat(inp, "d1", &d1)) d1 = 1. ;
    if (!sf_histfloat(inp, "o1", &o1)) o1 = 0. ;

    n2=sf_leftsize(inp,1);

    if (NULL != sf_getstring("zpo")) {
	zpo = sf_input("zpo");
	pos = sf_intalloc(n2);	
    } else {
	zpo = NULL;
	pos = NULL;
	sf_error("need the position file");
    }

    if (NULL != sf_getstring("zeropo")) {
	zeropo = sf_input("zeropo");
	zepo = sf_intalloc(n2);	
    } else {
	zeropo = NULL;
	zepo = NULL;
    }
    
    oeq=sf_floatalloc(n1);
    aeq=sf_floatalloc(n1);

    if(NULL != zpo) {
	sf_intread(pos,n2,zpo);
    }
    if(NULL != zeropo) {
	sf_intread(zepo,n2,zeropo);
    }
    
    for (i2=0; i2<n2; i2++) {
	sf_floatread(oeq,n1,inp);
	
	spo = 0;
	mspo = 0;
	
	sf_warning("slice %d of %d;",i2+1,n2);
	
	spo = pos[i2];
	mspo = spo+1;
	
	if(NULL != zeropo) {
	    zp = zepo[i2];
	    
	    for (i1=0; i1<n1; i1++) {
		if (i1 <= zp) {
		    aeq[i1] = 0.;
		} else if(i1 <= spo && i1 > zp) {
		    aeq[i1] = oeq[mspo];
		} else {
		    aeq[i1] = oeq[i1];
		}
	    }
	} else {
	    for (i1=0; i1<n1; i1++) {
		
		if(i1 <= spo) {
		    
		    aeq[i1] = oeq[mspo];
		} else {
		    aeq[i1] = oeq[i1];
		}
	    }	    
	}
	
	sf_floatwrite(aeq,n1,out);
    }
    
    exit(0);
}
