/* Bi-variate L1 regression */
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

#include "l1.h"

int main(int argc, char* argv[])
{
    int nd, n1, niter;
    float *d, *a, *b, alfa, beta, perc, fact;
    char *type;
    sf_file inp, reg, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    reg = sf_input("reg");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(reg,"n1",&n1) || n1 != nd)
	sf_error("Need n1=%d in reg",nd);
    if (!sf_histint(reg,"n2",&n1) || n1 != 2)
	sf_error("Need n2=2 in reg");
    if (NULL == (type = sf_getstring("type"))) type="threshold";
    /* thresholding type */

    sf_putint(out,"n1",2);
    
    d = sf_floatalloc(nd);
    a = sf_floatalloc(nd);
    b = sf_floatalloc(nd);

    sf_floatread(d,nd,inp);
    sf_floatread(a,nd,reg);
    sf_floatread(b,nd,reg);
    sf_fileclose(reg);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of POCS iterations */
    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening */
    if (!sf_getfloat("fact",&fact)) fact=1.5;
    /* factor for sharpening */

    l1_init(nd,niter,perc,fact,type,true);

    bil1(d,a,b,&alfa,&beta);
    
    sf_floatwrite(&alfa,1,out);
    sf_floatwrite(&beta,1,out);

    exit(0);
}
