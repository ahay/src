/* Convolution of two helix filters. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "conv.h"

int main(int argc, char* argv[]) 
{
    bool one;
    int na, nb, ns;
    sf_filter ss, aa, bb;
    char* file;
    sf_file inp, oth, out, lagin, lagout;

    sf_init(argc,argv);
    inp = sf_input("in");
    oth = sf_input("other");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&na)) sf_error("No n1= in input");
    aa = sf_allocatehelix(na);
    if (NULL == (file = sf_histstring(inp,"lag"))) sf_error("No lag= in input");
    lagin = sf_input(file);
    sf_intread(aa->lag,na,lagin);
    sf_floatread(aa->flt,na,inp);
    sf_fileclose(lagin);
    free(file);

    if (!sf_histint(oth,"n1",&nb)) sf_error("No n1= in other");
    bb = sf_allocatehelix(nb);
    if (NULL == (file = sf_histstring(oth,"lag"))) sf_error("No lag= in other");
    lagin = sf_input(file);
    sf_intread(bb->lag,nb,lagin);
    sf_floatread(bb->flt,nb,oth);
    free(file);

    if (!sf_getbool("one",&one)) one=true;
    /* include leading one */

    ss = conv (aa, bb, one);
    ns = ss->nh;

    if (NULL == (file = sf_getstring("lag"))) sf_error("Need lag=");
    lagout = sf_output(file);
    sf_settype(lagout,SF_INT);
    sf_putint(lagout,"n1",ns);
    sf_fileflush(lagout,lagin);
    sf_intwrite(ss->lag,ns,lagout);

    sf_putint(out,"n1",ns);
    sf_putstring(out,"lag",file);
    sf_floatwrite(ss->flt,ns,out);

    exit(0);
}
