/* Calculating frequency attenuation gradient. */
/*
  Copyright (C) 2011 Jilin University
  
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
#include "linearfit.h"

int main(int argc, char* argv[])
{
    bool sign, grad;
    int iw, nw, i1, n1, i2, n2, n1w, nd=0, wbeg=0, wend=0;
    float lperc, hperc, emax, freq;
    float d1, dw, w0, etotal, f1=0., f2=0., e1=0., e2=0., ecum, *fdg, *e;
/*    float cum1, cum2, cum3, cum4, numer, denom; */
    float *sol=NULL, *dat=NULL, **func=NULL;
    char *type;
    sf_complex *inp=NULL;
    sf_file in, out;
   
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    
    n2 = sf_leftsize(in,2);
    if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
    sf_unshiftdim(in, out, 2);
    n1w = n1*nw;

    if (!sf_getbool("grad",&grad)) grad=true;
    /* If y, output attenuation gradient; if n, output absorption factor */

    if (!sf_getfloat("lperc",&lperc)) lperc=65.;
    /* Low percentage of total energy */ 

    if (!sf_getfloat("hperc",&hperc)) hperc=85.;
    /* High percentage of total energy */ 

    if (lperc >= hperc || lperc < 0. || hperc > 100.) 
	sf_error("Need 0 <= lperc <= hperc <= 100.");

    if (NULL == (type=sf_getstring("type"))) type="attenuation";
    /* [low,full,ratio,attenuation] attribute type, the default is attenuation  */
    if(type[0]=='r') {
	if (!sf_getfloat("freq",&freq)) sf_error("Need freq when type=ratio");
        /* Frequency corresponding to energy ratio, valid when type=ratio */ 
    }
    
    if (SF_COMPLEX == sf_gettype(in)) {
	inp = sf_complexalloc(n1w);
	sf_settype(out,SF_FLOAT);
    }

    fdg = sf_floatalloc(n1);
    e   = sf_floatalloc(n1w);
 
    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);

	for (i1=0; i1 < n1; i1++) {
	    fdg[i1] = 0.;
	}
	if (SF_COMPLEX != sf_gettype(in)) {
	    sf_floatread(e,n1w,in);
	} else {
	    sf_complexread(inp,n1w,in);
	    for (i1=0; i1 < n1; i1++) {
		for (iw=0; iw < nw; iw++) {
		    e[iw*n1+i1] = sqrtf(crealf(inp[iw*n1+i1])*
					crealf(inp[iw*n1+i1])+
					cimagf(inp[iw*n1+i1])*
					cimagf(inp[iw*n1+i1]));
		}
	    }
	}
	for (i1=0; i1 < n1; i1++) {
	    switch (type[0]) {
		case 'a':
		    if (grad) {
			etotal = 0.;
			for (iw=0; iw < nw; iw++) {
			    etotal += e[iw*n1+i1];
			}
			ecum = 0.;
			sign = false;
			for (iw=0; iw < nw; iw++) {
			    ecum += e[iw*n1+i1];
			    if (!sign && ecum >= lperc*etotal/100.) {
				f1 = w0+iw*dw;
				e1 = e[iw*n1+i1];
				sign = true;
			    }
			    if (ecum >= hperc*etotal/100.) {
				f2 = w0+iw*dw;
				e2 = e[iw*n1+i1];
				break;
			    }
			}
			if (f1 == f2) {
			    fdg[i1] = 0.;
			} else {
			    fdg[i1] = (e2-e1)/(f2-f1);
			}
		    } else {
/*			numer = 0.;
			denom = 0.;
			
			for (iw=0; iw < nw; iw++) {
			    numer += e[iw*n1+i1]*(w0+iw*dw);
			    denom += e[iw*n1+i1];
			}
			
			index = (int)((numer/denom-w0)/dw);
			etotal = 0.;
			for (iw=0; iw < nw; iw++) {
			    etotal += e[iw*n1+i1];
			}
			ecum = 0.;
			for (iw=0; iw < nw; iw++) {
			    ecum += e[iw*n1+i1];
			    if (ecum >= hperc*etotal/100.) {
				wend = iw;
				break;
			    }
			}
			
			cum1 = 0.;
			cum2 = 0.;
			cum3 = 0.;
			cum4 = 0.;
			for (iw=index; iw < wend; iw++) {
			    cum1 += w0+iw*dw;
			    cum2 += logf(e[iw*n1+i1]);
			    cum3 += (w0+iw*dw)*logf(e[iw*n1+i1]);
			    cum4 += (w0+iw*dw)*(w0+iw*dw);
			}
			fdg[i1] = ((nw-index)*cum3-cum1*cum2)/
			    ((nw-index)*cum4-cum1*cum1);
*/
/*                      index = 0;                */
			emax = e[i1];
			for (iw=0; iw < nw-1; iw++) {
			    if (emax < e[(iw+1)*n1+i1]) {
				emax = e[(iw+1)*n1+i1];
				wbeg = iw+1;
			    }
			}
			
			wend = nw;
			
			if (wbeg == wend-1) wbeg = 0;
			
			nd = wend-wbeg+1;
			sol = sf_floatalloc(2);		
			dat = sf_floatalloc(nd);
			func = sf_floatalloc2(nd,2);
			
			for (iw=wbeg; iw < wend; iw++) {
			    dat[iw-wbeg] = logf(e[iw*n1+i1]);
			    func[0][iw-wbeg] = 1.;
			    func[1][iw-wbeg] = w0+iw*dw;
			}
			
			linearfit (nd,dat,func,sol);
			
			fdg[i1] = sol[1];
			
		    }
		    break;
		case 'l':
		    etotal = 0.;
		    for (iw=0; iw < nw; iw++) {
			etotal += e[iw*n1+i1];
		    }
		    ecum = 0.;
		    for (iw=0; iw < nw; iw++) {
			ecum += e[iw*n1+i1];
			if (ecum >= lperc*etotal/100.) {
			    f1 = w0+iw*dw;
			    break;
			}
		    }
		    fdg[i1] = f1;
		    break;
		case 'f':
		    etotal = 0.;
		    for (iw=0; iw < nw; iw++) {
			etotal += e[iw*n1+i1];
		    }
		    ecum = 0.;
		    for (iw=0; iw < nw; iw++) {
			ecum += e[iw*n1+i1];
			if (ecum >= hperc*etotal/100.) {
			    f2 = w0+iw*dw;
			    break;
			}
		    }
		    fdg[i1] = f2;
		    break;
		case 'r':
		    etotal = 0.;
		    for (iw=0; iw < nw; iw++) {
			etotal += e[iw*n1+i1];
		    }
		    ecum = 0.;
		    for (iw=0; iw < nw; iw++) {
			if(w0+iw*dw<=freq) {
			    ecum += e[iw*n1+i1];
			} else {
			    break;
			}
		    }
		    fdg[i1] = ecum/(etotal+FLT_EPSILON);
		    break;
		default:
		    sf_error("Unknown operator \"%s\"",type);
	    }
	}
	sf_floatwrite(fdg,n1,out);
    }
    sf_warning(".");
    
    exit(0);
}

/* 	$Id$	 */
