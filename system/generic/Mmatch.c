/* Simple matching filtering */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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
    bool adj;
    int n1, n2, i1, i2, i, j, nf;
    float *data, *noiz, *filt;
    sf_file inp, out, oth;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    oth = sf_input("other");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	/* input data, output filter */
	if (!sf_histint(inp,"n1",&n1)) sf_error("No n1=");
	if (!sf_histint(inp,"n2",&n2)) sf_error("No n2=");
	if (!sf_getint("nf",&nf)) sf_error("Need nf=");
	/* filter size */

	sf_putint(out,"n1",nf);
	sf_putint(out,"n2",1);
    } else {
	/* input filter, output data */
	if (!sf_histint(inp,"n1",&nf)) sf_error("No n1=");
	if (!sf_histint(oth,"n1",&n1)) sf_error("No n1=");
	if (!sf_histint(oth,"n2",&n2)) sf_error("No n2=");

	sf_fileflush(out,oth);  /* copy data dimensions */
    }

    filt = sf_floatalloc(nf);
    data = sf_floatalloc(n1);
    noiz = sf_floatalloc(n1);

    if (adj) {
	for (i=0; i < nf; i++) 
	    filt[i]=0;
    } else {
	sf_floatread(filt,nf,inp);
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(noiz,n1,oth);

	if (adj) {
	    sf_floatread(data,n1,inp); 
	} else {
	    for (i1=0; i1 < n1; i1++)
		data[i1] = 0.;
	}

	for (i=0; i < nf; i++) {
	    for (i1=0; i1 < n1; i1++) {
		j=i1-i+nf/2; /* symmetric filter */

		/* zero value boundary conditions */
		if (j < 0 || j >= n1) continue; 
		    
		if (adj) {
		    filt[i] += noiz[j]*data[i1]; 
		} else {
		    data[i1] += noiz[j]*filt[i];
		}
	    }
	}

	if (!adj) sf_floatwrite(data,n1,out);
    }

    if (adj) sf_floatwrite(filt,nf,out); 
		 
    exit(0);
}
