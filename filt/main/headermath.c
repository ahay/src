/* Mathematical operations, possibly on header keys. 

Known functions: cos,  sin,  tan,  acos,  asin,  atan, 
                 cosh, sinh, tanh, acosh, asinh, atanh,
                 exp,  log,  sqrt, abs

See also sfmath. 

An addition operation can be performed by sfstack.
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

int main(int argc, char* argv[])
{
    int i, i1, i2, n1, n2, len, nkey;
    sf_file in, out;
    char *eq, *output, *key, *arg;
    float **ftra, **fbuf, **fst, d2, o2;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    sf_putint(out,"N",0);
    sf_putint(out,"T",1);
    sf_putint(out,"input",2);

    for (i=0; i < SF_NKEYS; i++) {
	sf_putint(out,sf_segykeyword(i),i+3);
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	arg = argv[i];
	eq =  strchr(arg,'=');
	if (NULL == eq) continue; /* not a parameter */
	if (0 == strncmp(arg,"output",6)) continue; /* not a key */
	
	len = (size_t) (eq-arg);
	key = sf_charalloc(len+1);
	memcpy(key,arg,len);
	key[len]='\0';

	if (sf_getint(key,&nkey))
	    sf_putint(out,key,nkey+3);
	free(key);
    }

    if (SF_FLOAT != sf_gettype (in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) n1=1;
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (n1 > 1) { 
	if (n2 > 1) { /* input: many keys */
	    sf_putint(out,"n1",1);
	} else { /* input: one key, arranged in n1 */
	    n2 = n1;
	    n1 = 1;
	}
    }

    if (!sf_histfloat(in,n1>1? "d2":"d1",&d2)) d2=1.;
    if (!sf_histfloat(in,n1>1? "o2":"o1",&o2)) o2=0.;
    
    if (NULL == (output = sf_getstring("output"))) sf_error("Need output=");
    /* Describes the output in a mathematical notation. */

    len = sf_math_parse (output,out,SF_FLOAT);
    
    ftra = sf_floatalloc2(n1,n2);
    fbuf = sf_floatalloc2(n2,n1+3);
    fst  = sf_floatalloc2(n2,len+3);
    
    sf_floatread(ftra[0],n1*n2,in);

    for (i2=0; i2 < n2; i2++) {
	fbuf[0][i2]=(float) i2;  /* N */
	fbuf[1][i2]=o2+i2*d2;    /* T */
	fbuf[2][i2]=ftra[0][i2]; /* input */ 
    }
    for (i1=0; i1 < n1; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    fbuf[i1+3][i2]=ftra[i2][i1];
	}
    }

    free (ftra[0]);
    free (ftra);

    sf_math_evaluate (len, n2, fbuf, fst);
    sf_floatwrite(fst[1],n2,out);

    exit(0);
}

/* 	$Id$	 */
