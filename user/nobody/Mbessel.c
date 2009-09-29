/* Bessel functions */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "bessel_I0.h"
#include "bessel_I1.h"

int main(int argc, char* argv[])
{
    off_t nsiz;
    int nbuf, k, order;
    char *type;
    float *fbuf;
    sf_file inp, out;
    double (*func)(double);

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (NULL == (type = sf_getstring("type"))) type="J";
    if (!sf_getint("order",&order)) order=0;

    switch(type[0]) {
	case 'i':
	case 'I':
	    switch (order) {
		case 0:
		    func = bessel_I0;
		    break;
		case 1:
		    func = bessel_I1; 
		    break;
		default:
		    sf_error("Order %d for type %s is not implemented",order,type);
	    }
	    break;
	default:
	    sf_error("Type %s is not implemented",type);
	    break;
    }

    nbuf = BUFSIZ/sizeof(float);
    fbuf = sf_floatalloc(nbuf);

    for (nsiz = sf_filesize(inp); nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;
	sf_floatread(fbuf,nbuf,inp);
	for (k=0; k < nbuf; k++) {
	    fbuf[k] = func(fbuf[k]);
	}
	sf_floatwrite(fbuf,nbuf,out);
    }

    exit(0);
}
