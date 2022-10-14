/* Estimate equivalent Q value based on a reference layer by using Local centroid frequency shift (LCFS) method. */
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

int main (int argc,char *argv[])
{
    bool verb;

    float o1, d1;
    int n1, n2, i1, i2, i;
    int spo, mspo, tp, *repo=NULL;
    float detf, dif;
    float *cf, *var=NULL, *eqt;
    sf_file inp, out, var2=NULL, repos=NULL;

    sf_init (argc,argv);

    inp = sf_input("in");
    out = sf_output("out");    

    if (!sf_getbool("verb",&verb)) verb=false;
	
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");

    n2=sf_leftsize(inp,1);

    if (NULL != sf_getstring("var2")) { /* variance */
	var2 = sf_input("var2");
	var = sf_floatalloc(n1);
    } else {
	sf_error("need var2.");
    }

    if (NULL != sf_getstring("repos")) { /* Position of reference point */
	repos = sf_input("repos");
	repo = sf_intalloc(n2);
    } else {
	sf_error("need repos.");
    }

    cf = sf_floatalloc(n1);    
    eqt = sf_floatalloc(n1);

    if(NULL != repos) {
	sf_intread(repo,n2,repos);
    }
    
    for (i2=0; i2<n2; i2++) {
	sf_floatread(cf,n1,inp);

	if(NULL != var2) {
	    sf_floatread(var,n1,var2);
	}
	
	spo = 0;
	spo = repo[i2];
	mspo = spo+1;
	
	for (i1=mspo; i1<n1; i1++) {
	    tp = 0;
	    tp = i1;
	    detf = 0.;
	    for (i=mspo; i<=tp; i++) {
		dif = cf[i-1]-cf[i];
		if (dif < 0) dif = 0.000001;
		detf += (dif)/var[i-1] ;
	    }	    
	    eqt[i1] = SF_PI*(i1-spo)*d1/detf;
	}

	for (i1=0; i1<mspo; i1++) {
	    eqt[i1] = eqt[mspo];
	}
	
	sf_floatwrite(eqt,n1,out);
    }
  
}
