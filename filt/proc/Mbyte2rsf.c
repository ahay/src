/* Convert raw byte images to RSF.

Takes: < raw.img > out.rsf
*/
#include <stdio.h>

#include <rsf.h> 

int main(int argc, char* argv[])
{ 
    int n1, n2, x, y;
    unsigned char *line;
    float* array;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");

    if(!sf_getint("n1",&n1)) sf_error("Need n1="); 
    /* vertical dimension */
    if(!sf_getint("n2",&n2)) sf_error("Need n2="); 
    /* horizontal dimension */

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",1.);
    sf_putfloat(out,"d2",1.);
    sf_setformat(out,"native_float");

    array = sf_floatalloc(n1);
    line = (unsigned char *) sf_alloc(n1,sizeof(unsigned char));

    for (y = 0; y < n2; y++) {
	if (n1 != fread(line, sizeof(unsigned char), n1, stdin))
	    sf_error("trouble reading input data");
	for (x = 0; x < n1; x++) {
	    array[x] = (float) line[x];
	}
	sf_floatwrite(array,n1,out); 
    }

    exit (0);
}

