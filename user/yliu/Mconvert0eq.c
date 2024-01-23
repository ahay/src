/* Convert equivalent Q value from reference layer to t0 location. 

Ignore this step if selecting t0 as the reference layer in module "sflcfs".
*/
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
    int   spo, mspo, *repo;
    float d1, o1;
    float *eq, *teq;

    sf_file inp, out, repos;

    sf_init(argc, argv);
    inp  = sf_input("in");
    out  = sf_output("out");
    sf_settype(out,SF_FLOAT);

    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;

    if (!sf_histint(inp, "n1", &n1)) sf_error("No n1=in input");
    if (!sf_histfloat(inp, "d1", &d1)) d1 = 1. ;
    if (!sf_histfloat(inp, "o1", &o1)) o1 = 0. ;

    n2=sf_leftsize(inp,1);
    
    eq = sf_floatalloc(n1);
    teq = sf_floatalloc(n1);

    if (NULL != sf_getstring("repos")) { /* Position of reference point */
	repos = sf_input("repos");
	repo = sf_intalloc(n2);
    } else {
	repos = NULL;
	repo = NULL;
	sf_error("need repos.");
    }

    if(NULL != repos) {
	sf_intread(repo,n2,repos);
    }

    for (i2=0; i2<n2; i2++) {
	
	sf_floatread(eq, n1, inp);
	/* read equivalent Q */

	spo = 0;
	spo = repo[i2];
	mspo = spo+1;
	
	for (i1 = 0; i1 < n1; i1++) {
	    teq[i1] = 0.;
	}
	for (i1 = 0; i1 < n1; i1++) {

	    if (i1 < mspo) {
		teq[i1] = eq[i1];
	    } else {
		teq[i1] = i1*d1/((i1-spo)*d1/eq[i1]+spo*d1/eq[spo]);
	    }
	}
	sf_floatwrite(teq,n1,out);
    }
}
