/* Extend a dataset by duplicating in the second axis dimension.

Takes: < input.rsf > output.rsf

This operation is adjojnt to sfstack.
*/

#include <string.h>
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int j, n1, n2, n3, i2, i3, ni, esize;
    float f;
    size_t n;
    sf_file in, out;
    char key1[7], key2[7], *val, *trace;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"esize",&esize)) sf_error("No esize= in input");

    n = (size_t) (n1*esize);

    if (!sf_getint("n2",&n2)) sf_error("Need n2=");
    /* Size of the newly created dimension */ 
    sf_putint(out,"n2",n2);
    if (sf_getfloat("d2",&f)) sf_putfloat(out,"d2",f);
    /* Sampling of the newly created dimension */ 
    if (sf_getfloat("o2",&f)) sf_putfloat(out,"o2",f);
    /* Origin of the newly created dimension */
    
    n3 = 1;
    for (j=2; j < SF_MAX_DIM; j++) {
	sprintf(key2,"n%d",j+1);
	sprintf(key1,"n%d",j);
	if (!sf_histint(in,key1,&ni)) break;
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j+1);
	sprintf(key1,"o%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j+1);
	sprintf(key1,"d%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j+1);
	sprintf(key1,"label%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }
    
    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    trace = sf_charalloc (n);
    
    for (i3=0; i3 < n3; i3++) {
	sf_read (trace, 1, n, in);
	for (i2=0; i2 < n2; i2++) {
	    sf_write(trace, 1, n, out);
	} 
    }
    
    exit (0);
}

/* 	$Id: spray.c,v 1.3 2003/09/29 14:34:56 fomels Exp $	 */
