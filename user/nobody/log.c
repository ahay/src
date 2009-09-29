/* Take the base-10 logarithm of input data.

Alternatively, use sfmath or sfheadermath.
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    float range, *data, di, avg=0., big;
    int esize, centered;
    size_t i, size;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float data type");
    if(!sf_histint(in,"esize",&esize)) esize=4;
    size = (size_t) sf_filesize (in);

    if(!sf_getint("centered",&centered)) centered=2;
    /* [0,1,2] defines method of shifting mean */
    if(!sf_getfloat("range",&range)) range=3.;
    /* Smallest allowed value */

    data = sf_floatalloc(size);
    sf_floatread (data,size,in);

    big = (float) FLT_MIN_10_EXP;
    for (i=0; i< size; i++) {
	di = fabsf(data[i]);
	di = (0. >= di)? (float) FLT_MIN_10_EXP: log10f(di);
	if (di > big) big=di;
	data[i] = di;
    }

    for (i=0; i< size; i++) {
	di = data[i];
	if (di < big-range) di = big-range;
    }

    switch (centered) {
	case 0:
	    avg=big;
	    sf_putstring(out,"label2","log base 10");
	    break;
	case 1:
	    sf_putint(out,"pclip",100);
	    avg = (big + (big-range)) / 2.;
	    break;
	case 2:
	    avg = 0.;
	    for (i=0; i< size; i++) {
		avg += data[i];
	    }
	    avg /= size;
	    sf_putfloat(out,"range",2.*(big-avg));
	    break;
	default:
	    sf_error ("Unsupported centered=%d",centered);
	    break;
    }

    for (i=0; i< size; i++) {
	data[i] -= avg;
    }

    sf_floatwrite (data,size,out);


    exit (0);
}
