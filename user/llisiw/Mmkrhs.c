/* Take out traveltime according to receiver location */
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
    bool elev;
    int n[SF_MAX_DIM], it, nt, dim, i, is, nshot, **m, rhs, nrhs;
    float **f, **e;
    sf_file in, mask, out;
    
    sf_init(argc,argv);
    in   = sf_input("in");
    out  = sf_output("out");
    mask = sf_input("mask");
    
    /* read input dimension */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	nt *= n[i];
    }

    /* read number of shots */
    if (!sf_histint(in,"n4",&nshot)) nshot=1;
    nt /= nshot;

    /* read complete timetable */
    f = sf_floatalloc2(nt,nshot);
    sf_floatread(f[0],nt*nshot,in);

    /* read receiver location */
    m = sf_intalloc2(nt,nshot);
    sf_intread(m[0],nt*nshot,mask);

    nrhs = 0;
    for (it=0; it < nt; it++) {
	if (m[0][it]==1) 
	    nrhs++;
    }

    e = sf_floatalloc2(nrhs,nshot);

    sf_putint(out,"n1",nrhs);
    sf_putint(out,"n2",nshot);
    sf_putint(out,"n3",1);
    sf_putint(out,"n4",1);

    /* extract right-hand side from complete timetable */
    for (is=0; is < nshot; is++) {
	
	rhs = 0;
	for (it=0; it < nt; it++) {
	    if (m[is][it]==1) {
		e[is][rhs] = f[is][it];
		rhs++;
	    }
	}
    }

    sf_floatwrite(e[0],nrhs*nshot,out);

    exit(0);
}
