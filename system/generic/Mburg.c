/* Burg's method for 1-D PEF estimation */
/*
  Copyright (C) 2022 The University of Texas at Austin

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

#include "burg.h"

int main(int argc, char* argv[]) 
{
    int nd, na, i, dim, n[SF_MAX_DIM];
    float *d, *dpef, *pef;
    char key[5];
    sf_file inp, out, filt;
    
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    filt = sf_output("filter");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    dim = sf_filedims(inp,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }

    if (!sf_getint("na",&na)) na=10;
    /* filter size */

    sf_putint(filt,"n1",na);
    for (i=1; i < dim; i++) {
	snprintf(key,3,"n%d",i+1);
	sf_putint(filt,key,1);
    }
    
    d = sf_floatalloc(nd);
    dpef = sf_floatalloc(nd);
    pef = sf_floatalloc(na);

    burg_init(nd,na);

    sf_floatread(d,nd,inp);
    
    burg(d,dpef,pef);

    sf_floatwrite(dpef,nd,out);
    sf_floatwrite(pef,na,filt);

    exit(0);
}
