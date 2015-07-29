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
    int i, i2, n1, nbuf, **buf, *buf1, nk;
    float o1, d1;
    off_t n2, nleft;
    const char *key;
    char *arg;
    sf_file in, keys[SF_MAXKEYS], out, tfile;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint(in,"n1",&n1) &&
	!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* number of samples in a trace */
    if (!sf_histfloat(in,"d1",&d1) &&
	!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* trace sampling */
    if (!sf_histfloat(in,"o1",&o1) &&
	!sf_getfloat("o1",&o1)) o1=0;
    /* trace origin */

    n2 = sf_leftsize(in,1);

    nbuf = BUFSIZ/sizeof(int);

    if (NULL != sf_getstring("tfile")) {
	tfile = sf_input("tfile"); /* trace header file */
	if (SF_INT != sf_gettype(tfile))
	    sf_error("Need integer data in tfile");
	if (!sf_histint(tfile,"n1",&nk) || (SF_NKEYS > nk))
	    sf_error ("Need at least n1=%d keys in tfile",SF_NKEYS);
	if (nk*n2 != sf_filesize(tfile))
	    sf_error ("Wrong number of traces in tfile");
    } else {
	tfile = NULL;
	nk = SF_NKEYS;
    }

    sf_putint(out,"n1",nk);
    sf_settype(out,SF_INT);

    if (NULL != tfile) sf_fileflush(out,tfile);

    buf = sf_intalloc2(nk,nbuf);
    buf1 = sf_intalloc(nbuf);

    segy_init(nk,tfile);

    for (i=0; i < nk; i++) {
	key = segykeyword(i);
	if (NULL != (arg = sf_getstring(key))) {
	    keys[i] = sf_input(key);
	    if (SF_INT != sf_gettype(keys[i]))
		sf_error("Need integer data in file \"%s\"",arg); 
	    if (n2 != sf_filesize(keys[i])) 
		sf_error("Need filesize=%lld in file \"%s\"",n2,arg); 
	    free(arg);
	} else {
	    keys[i] = NULL;
	    for (i2=0; i2 < nbuf; i2++) {
		buf[i2][i] = 0;
	    }
	}
    }

    for (nleft=n2; nleft > 0; nleft -= nbuf) {
	if (nbuf > nleft) nbuf = nleft;

	/* read from initial trace header file */
	if (NULL != tfile) sf_intread(buf[0],nk*nbuf,tfile);

	for (i=0; i < nk; i++) {
	    key = segykeyword(i);

	    if (NULL != keys[i]) {
		sf_intread(buf1,nbuf,keys[i]);
		for (i2=0; i2 < nbuf; i2++) {
		    buf[i2][i] = buf1[i2];
		}
	    } else { /* change ns, dt, and delrt */
		if (0==strcmp(key,"ns")) {
		    for (i2=0; i2 < nbuf; i2++) {
			buf[i2][i] = n1;
		    }
		} else if (0==strcmp(key,"dt")) {
		    for (i2=0; i2 < nbuf; i2++) {
			buf[i2][i] = (int) (d1*1000000. + 0.5);
		    }
		} else if (0==strcmp(key,"delrt") && o1 != 0) {
		    keys[i] = NULL;
		    for (i2=0; i2 < nbuf; i2++) {
			buf[i2][i] = (o1>0)? (int) (o1*1000. + 0.5): (int) (o1*1000. - 0.5);
		    }
		}
	    }
	}
	
	sf_intwrite(buf[0],nk*nbuf,out);
    }

    free(buf1);
    free(*buf); free(buf);

    exit(0);
}
