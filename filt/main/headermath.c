/* Mathematical operations, possibly on header keys.

Takes: < input.rsf > output.rsf

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

    for (i=0; i < SF_NKEYS; i++) {
	sf_putint(out,sf_segykeyword(i),i+2);
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
	    sf_putint(out,key,nkey+2);
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

    len = sf_math_parse (output,out);
    
    ftra = sf_floatalloc2(n1,n2);
    fbuf = sf_floatalloc2(n2,n1+2);
    fst  = sf_floatalloc2(n2,len+2);
    
    sf_read(ftra[0],sizeof(float),n1*n2,in);

    for (i2=0; i2 < n2; i2++) {
	fbuf[0][i2]=(float) i2;
	fbuf[1][i2]=o2+i2*d2;
    }
    for (i1=0; i1 < n1; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    fbuf[i1+2][i2]=ftra[i2][i1];
	}
    }

    free (ftra[0]);
    free (ftra);

    sf_math_evaluate (len, n2, fbuf, fst);
    sf_write(fst[1],sizeof(float),n2,out);
    
    exit(0);
}

/* 	$Id: headermath.c,v 1.5 2003/09/29 14:34:56 fomels Exp $	 */

