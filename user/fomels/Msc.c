/* Surface-consistent decomposition */
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

#include "sc.h"

int main(int argc, char* argv[])
{
    bool prec;
    int nd, nm, ns, nx, n2, im, min, max, id, i, niter;
    int **indx, *size;
    float *model, *data;
    char name[6];
    sf_file inp, index, out;

    sf_init(argc,argv);
    
    inp = sf_input("in");
    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");

    index = sf_input("index");
    if (SF_INT != sf_gettype(index)) sf_error("Need integet index");

    if (!sf_histint(index,"n1",&nm)) sf_error("No n1= in index");
    if (!sf_histint(index,"n2",&n2) || n2 != nd) sf_error("Need n2=%d in index",nd);

    data = sf_floatalloc(nd);
    indx = sf_intalloc2(nm,nd);
    size = sf_intalloc(nm);

    if (!sf_getint("niter",&niter)) niter=0; /* number of iterations */
    if (!sf_getbool("prec",&prec)) prec=true; /* if apply preconditioning */    

    sf_intread(indx[0],nd*nm,index);

    nx = 0;
    for (im=0; im < nm; im++) {
	min = max = indx[0][im];
	for (id=1; id < nd; id++) {
	    i = indx[id][im];
	    if (i < min) min=i;
	    if (i > max) max=i;
	}
	size[im]=min-max;
	nx += size[im];
    }
    model = sf_floatalloc(nx);

    sc_init(nm, indx, size);
    
    sf_floatread(data,nd,inp);

    sf_solver(sc_lop,sf_cgstep,nx,nd,model,data,niter,"verb",true,"end");

    nx = 0;
    for (im=0; im < nm; im++) {
	snprintf(name,6,"out%d",im+1);
	out = sf_output(name);
	ns = size[im];
	sf_putint(out,"n1",ns);
	sf_floatwrite(model+nx,ns,out);
	sf_fileclose(out);
	nx += size[im];
    }
    
    exit(0);
}
    
    
