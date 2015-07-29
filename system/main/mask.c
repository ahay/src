/* Create a mask.

Mask is an integer data with ones and zeros. 
Ones correspond to input values between min and max.

The output can be used with sfheaderwindow.
*/
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
#include <float.h>
#include <limits.h>

#include <rsf.h>

int main(int argc, char* argv[]) {
    int *ibuf, imin=-INT_MAX, imax=INT_MAX;
    size_t nbuf, j, bufsiz;
    off_t nsiz;
    float *fbuf=NULL, fmin=-FLT_MAX, fmax=FLT_MAX;
    sf_datatype type;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    type = sf_gettype(in);
    bufsiz = sf_bufsiz(in);
    sf_settype(out,SF_INT);

    switch(type) {
	case SF_FLOAT:
	    nbuf = bufsiz/sizeof(float);
	    fbuf = sf_floatalloc (nbuf);

	    sf_getfloat("min",&fmin);
	    /* minimum header value */
	    sf_getfloat("max",&fmax);
	    /* maximum header value */
	    break;
	case SF_INT:
	    nbuf = bufsiz/sizeof(int);

	    sf_getint("min",&imin);
	    /* minimum header value */
	    sf_getint("max",&imax);
	    /* maximum header value */
	    break;
	default:
	    nbuf = 0;
	    sf_error("Unsupported type %d",type);
	    break;
    }

    ibuf = sf_intalloc (nbuf);

    for (nsiz = sf_filesize (in); nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;

	switch(type) {
	    case SF_FLOAT:
		sf_floatread(fbuf,nbuf,in);
		for (j=0; j < nbuf; j++) {
		    ibuf[j] = (fbuf[j] <= fmax && fbuf[j] >= fmin);
		}
		break;
	    case SF_INT:
		sf_intread(ibuf,nbuf,in);
		for (j=0; j < nbuf; j++) {
		    ibuf[j] = (ibuf[j] <= imax && ibuf[j] >= imin);
		}
		break;
	    default:
		sf_error("Unsupported type %d",type);
		break;
	}

	sf_intwrite(ibuf,nbuf,out);
    }


    exit(0);
}
	    
/* 	$Id: mask.c 7791 2011-10-29 13:22:26Z sfomel $	 */
