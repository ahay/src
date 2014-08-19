/* Mathematical operations, possibly on header keys. 

Known functions for float data: 
cos,  sin,  tan,  acos,  asin,  atan, 
cosh, sinh, tanh, acosh, asinh, atanh,
exp,  log,  sqrt, abs, erf, erfc, sign

Known functions for int data: sign, abs

See also sfmath.

An addition operation can be performed by sfadd.
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
    bool segy;
    int i, i1, i2, n1, n2, n3, n, nt, len, nkey;
    sf_file in, out;
    int mem; /* for avoiding int to off_t typecast warning */
    off_t memsize;
    char *eq=NULL, *output=NULL, *key=NULL, *arg=NULL;
    float **ftra=NULL, **fbuf=NULL, **fst=NULL, d2, o2;
    int **itra=NULL, **ibuf=NULL, **ist=NULL;
    sf_datatype type;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    sf_putint(out,"N",0);
    sf_putint(out,"T",1);
    sf_putint(out,"input",2);

    type = sf_gettype(in);

    if (SF_FLOAT != type && SF_INT != type) sf_error("Need float or int input");

    if (!sf_histint(in,"n1",&n1)) n1=1;
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2); /* left dimensions after the first two */
    if (n1 > 1) {
	if (n2 > 1) { /* input: many keys */
	    sf_putint(out,"n1",1);
	} else { /* input: one key, arranged in n1 */
	    n2 = n1;
	    n1 = 1;
	}
    }

    if (!sf_getbool("segy",&segy)) segy=true;
    /* if SEGY headers */

    if (segy) {
	segy_init(n1,in);
    } else {
	other_init(n1,in);
    }

    for (i=0; i < n1; i++) {
	sf_putint(out,segykeyword(i),i+3);
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	arg = argv[i];
	eq =  strchr(arg,'=');
	if (NULL == eq) continue; /* not a parameter */
	if (0 == strncmp(arg,"output",6) ||
	    0 == strncmp(arg,    "--",2)) continue; /* not a key */

	len = (size_t) (eq-arg);
	key = sf_charalloc(len+1);
	memcpy(key,arg,len);
	key[len]='\0';

	if (sf_getint(key,&nkey))
	    sf_putint(out,key,nkey+3);
	free(key);
    }



    if (!sf_histfloat(in,n1>1? "d2":"d1",&d2)) d2=1.;
    if (!sf_histfloat(in,n1>1? "o2":"o1",&o2)) o2=0.;

    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");
    /* Describes the output in a mathematical notation. */

    if (!sf_getint("memsize",&mem))
        mem=sf_memsize();
    /* Max amount of RAM (in Mb) to be used */
    memsize = mem * (1<<20); /* convert Mb to bytes */

    len = sf_math_parse (output,out,type);

    /* number of traces for optimal I/O */
    nt = SF_MAX(1,memsize/((2*n1+len+6)*sizeof(float)));

    if (SF_FLOAT == type) { /* float type */
	ftra = sf_floatalloc2(n1,nt);
	fbuf = sf_floatalloc2(nt,n1+3);
	fst  = sf_floatalloc2(nt,len+3);
    } else {               /* int type */
	itra = sf_intalloc2(n1,nt);
	ibuf = sf_intalloc2(nt,n1+3);
	ist  = sf_intalloc2(nt,len+3);
    }

    for (n=n2*n3; n > 0; n -= nt) {
	if (n < nt) nt=n;

	if (SF_FLOAT == type) { 
	    sf_floatread(ftra[0],n1*nt,in);
	} else {
	    sf_intread(itra[0],n1*nt,in);
	}

	for (i2=0; i2 < nt; i2++) {
	    if (SF_FLOAT == type) { 
		fbuf[0][i2]=(float) i2;  /* N */
		fbuf[1][i2]=o2+i2*d2;    /* T */
		fbuf[2][i2]=ftra[0][i2]; /* input */
	    } else {
		ibuf[0][i2]=i2;          /* N */
		ibuf[1][i2]=o2+i2*d2;    /* T */
		ibuf[2][i2]=itra[0][i2]; /* input */
	    }
	}
	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2 < nt; i2++) {
		if (SF_FLOAT == type) { 
		    fbuf[i1+3][i2]=ftra[i2][i1];
		} else {
		    ibuf[i1+3][i2]=itra[i2][i1];
		}
	    }
	}
	
	if (SF_FLOAT == type) { 
	    sf_math_evaluate (len, nt, fbuf, fst);
	    sf_floatwrite(fst[1],nt,out);
	} else {
	    sf_int_math_evaluate (len, nt, ibuf, ist);
	    sf_intwrite(ist[1],nt,out);
	}
    }

    exit(0);
}
