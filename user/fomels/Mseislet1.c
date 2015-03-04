/* 1-D seislet transform */
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

#include <rsf.h>

#include "seislet1.h"

int main(int argc, char *argv[])
{
    int n1, n2, i2;
    bool inv, adj, unit, verb;
    char *type;
    sf_complex *pp, *qq; 
    float *dd;
    sf_file in, out, dip;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("freq"); /* variable frequency */

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    pp = sf_complexalloc(n1);
    qq = sf_complexalloc(n1);
    dd = sf_floatalloc(n1);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (NULL == (type=sf_getstring("type"))) type="haar";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    seislet1_init(n1,inv,unit,type[0]);
    seislet1_set(dd);

    for (i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("slice %d of %d;",i2+1,n2);

	sf_complexread(pp,n1,in);
	sf_floatread(dd,n1,dip);

	if (adj) {
	    seislet1_lop(adj,false,n1,n1,qq,pp);
	} else {
	    seislet1_lop(adj,false,n1,n1,pp,qq);
	}
	sf_complexwrite(qq,n1,out);
    }
    sf_warning(".");
    exit(0);
}
