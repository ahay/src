/* Make a trace header file for segywrite.

   Use the output for tfile= argument in segywrite.
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

#include <string.h>
#include <rsf.h>
#include "segy.h"

int main(int argc, char* argv[])
{
    int i, i2, n1, nbuf, *buf[SF_NKEYS], buf2[SF_NKEYS];
    float d1;
    sf_file in=NULL, keys[SF_NKEYS], out=NULL;
    off_t n2, nleft;
    char *key=NULL, *arg=NULL, zero[BUFSIZ];

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint(in,"n1",&n1) &&
	!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* number of samples in a trace */
    if (!sf_histfloat(in,"d1",&d1) &&
	!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* trace sampling */

    n2 = sf_leftsize(in,1);
    sf_putint(out,"n1",SF_NKEYS);
    sf_settype(out,SF_INT);

    nbuf = BUFSIZ/sizeof(int);
    memset(zero,0,BUFSIZ);

    for (i=0; i < SF_NKEYS; i++) {
	key = segykeyword(i);
	if (0==strcmp(key,"ns")) {
	    keys[i] = NULL;
	    buf[i] = sf_intalloc(nbuf);
	    for (i2=0; i2 < nbuf; i2++) {
		buf[i][i2] = n1;
	    }
	} else if (0==strcmp(key,"dt")) {
	    keys[i] = NULL;
	    buf[i] = sf_intalloc(nbuf);
	    for (i2=0; i2 < nbuf; i2++) {
		buf[i][i2] = (int) (d1*1000000.);
	    }
	} else if (NULL != (arg = sf_getstring(key))) {
	    keys[i] = sf_input(key);
	    if (SF_INT != sf_gettype(keys[i]))
		sf_error("Need integer data in file \"%s\"",arg); 
	    if (n2 != sf_filesize(keys[i])) 
		sf_error("Need filesize=%lld in file \"%s\"",n2); 
	    free(arg);

	    buf[i] = sf_intalloc(nbuf);
	} else {
	    keys[i] = NULL;
	    buf[i] = (int*) zero;
	}
    }

    for (nleft=n2; nleft > 0; nleft -= nbuf) {
	if (nbuf > nleft) nbuf = nleft;

	
	for (i=0; i < SF_NKEYS; i++) {
	    if (NULL != keys[i]) {
		sf_intread(buf[i],nbuf,keys[i]);
	    }
	}

	for (i2=0; i2 < nbuf; i2++) {
	    for (i=0; i < SF_NKEYS; i++) {
		buf2[i] = buf[i][i2];
	    }
	    sf_intwrite(buf2,SF_NKEYS,out);
	}
    }

    exit(0);
}
