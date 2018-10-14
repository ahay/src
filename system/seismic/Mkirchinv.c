/* Kirchhoff 2-D post-stack least-squares time migration with antialiasing. 

 Antialiasing by reparameterization. */
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

#include "kirchnew.h"
#include "tfweight.h"
#include "invert.h"

int main(int argc, char* argv[])
{
    int n12, n1, n2, n3, i3, sw, nw, niter;
    bool hd, ps;
    float *data, *modl, *modl0, *vrms, *error=NULL, *ww, *ff;
    float o1,d1,o2,d2;
    char *errfile;
    sf_file in, out, vel, in0, err=NULL, fwght, wght;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (NULL != sf_getstring("fweight")) {
	fwght = sf_input("fweight");
	nw = kiss_fft_next_fast_size((n1+1)/2)+1;
	ff = sf_floatalloc(nw*n2);
	sf_floatread(ff,nw*n2,fwght);
	sf_fileclose(fwght);

	wght = sf_input("weight");
	ww = sf_floatalloc(n1*n2);
	sf_floatread(ww,n1*n2,wght);
	sf_fileclose(wght);
	   
	tfweight_init(n1,nw,n2,ww,ff);
    } else {
	fwght = NULL;
	wght = NULL;
    }
    

    if (!sf_getbool("hd",&hd)) hd=true;
    /* if y, apply half-derivative filter */
    if (!sf_getbool("ps",&ps)) ps=true;
    /* if y, apply pseudo-unitary weighting */
    if (!sf_getint("sw",&sw)) sw=0;
    /* if > 0, select a branch of the antialiasing operation */
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    vrms = sf_floatalloc(n1);

    vel = sf_input("velocity");
    sf_floatread(vrms,n1,vel);
    sf_fileclose(vel);
    
    n12 = n1*n2;
    data = sf_floatalloc(n12);
    modl = sf_floatalloc(n12);

    if (NULL != sf_getstring("model0")) {
        in0 = sf_input("model0");
        modl0 = sf_floatalloc(n12);
        sf_floatread(modl0,n12,in0);
        sf_fileclose(in0);
    } else {
        modl0 = NULL;
    } 

    if (niter > 0) {
	errfile = sf_getstring("err");
	/* output file for error */
	if (NULL != errfile) {
	    err = sf_output(errfile);
	    sf_putint(err,"n1",niter);
	    sf_putfloat(err,"d1",1);
	    sf_putfloat(err,"o1",1);
	    sf_putstring(err,"label1","Iteration Number");
	    sf_putstring(err,"label2","Relative Squared Error");
	    sf_putint(err,"n2",1);
	}
	error = sf_floatalloc(niter);
    }

    kirchnew_init (vrms, o1, d1, d2, n1, n2, sw, ps, hd);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread (data,n12,in);

	if (NULL != fwght) {
	    invert(kirchnew_lop,tfweight_lop,
		   niter,n12,n12,modl,modl0,data,error);
	} else {
	    invert(kirchnew_lop,NULL,niter,n12,n12,modl,modl0,data,error);
	}

	sf_floatwrite (modl,n12,out);
	if (NULL != err) sf_floatwrite(error,niter,err);
    }


    exit(0);
}
