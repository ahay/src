/* Make receiver mask according to elevation info */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
    int n[SF_MAX_DIM], n1, n2, nt, dim, **m, nshot, ns, i, it, temp, index;
    float o[SF_MAX_DIM], o1, o2, d[SF_MAX_DIM], d1, d2, **ele;
    char key[4];
    sf_file in, elev, out;

    sf_init(argc,argv);
    in   = sf_input("in");
    out  = sf_output("out");
    elev = sf_input("elev");

    /* read input dimension */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }

    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
	sf_putint(out,"n3",1);
	sf_putfloat(out,"o3",o[1]);
	sf_putfloat(out,"d3",d[1]);
    }

    /* read elevation file dimentions */
    if (!sf_histint(elev,"n1",&n1)) n1 = n[1];
    if (!sf_histint(elev,"n2",&n2)) n2 = n[2];

    if (!sf_histfloat(elev,"o1",&o1)) o1 = o[1];
    if (!sf_histfloat(elev,"o2",&o2)) o2 = o[2];

    if (!sf_histfloat(elev,"d1",&d1)) d1 = d[1];
    if (!sf_histfloat(elev,"d2",&d2)) d2 = d[2];

    /* read number of shots */
    if (!sf_histint(elev,"n3",&nshot)) nshot=1;

    /* write output dimentions */
    sf_putint(out,"n4",nshot);
    sf_putfloat(out,"o4",1.);
    sf_putfloat(out,"d4",1.);
    sf_settype(out,SF_INT);

    m = sf_intalloc2(nt,nshot);
    ele = sf_floatalloc2(n1*n2,nshot);

    sf_floatread(ele[0],n1*n2*nshot,elev);

    for (ns=0; ns < nshot; ns++) {
	for (it=0; it < nt; it++) {
	    m[ns][it] = 0;
	}

	for (it=0; it < n1*n2; it++) {
	    if (ele[ns][it] < o[0]) 
		continue;
	    else {
		temp = it/n1;
		index = (temp*d2+o2-o[2])/d[2]*n[0]*n[1]
		    +((it-temp*n1)*d1+o1-o[1])/d[1]*n[0]
		    +(ele[ns][it]-o[0])/d[0];
		m[ns][index] = 1;
	    }
	}
    }

    sf_intwrite(m[0],nt*nshot,out);

    exit(0);
}
