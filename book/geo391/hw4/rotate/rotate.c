/* Rotate around. */
#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2, k1, k2;
    float x1, x2, c1, c2, cosa, sina, angle;
    float **orig, **rotd;
    char *interp;
    sf_file inp, out;

    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in inp");

    /* get parameters from command line */
    if (!sf_getfloat("angle",&angle)) angle=90.;
    /* rotation angle */

    if (NULL == (interp = sf_getstring("interp"))) 
	interp="nearest";
    /* [n,l,c] interpolation type */

    /* convert degrees to radians */
    angle *= SF_PI/180.;
    cosa = cosf(angle);
    sina = sinf(angle);

    orig = sf_floatalloc2(n1,n2);
    rotd = sf_floatalloc2(n1,n2);

    /* read data */
    sf_floatread(orig[0],n1*n2,inp);

    /* central point */
    c1 = (n1-1)*0.5;
    c2 = (n2-1)*0.5;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {

	    /* rotated coordinates */
	    x1 = c1+(i1-c1)*cosa-(i2-c2)*sina;
	    x2 = c2+(i1-c1)*sina+(i2-c2)*cosa;

	    /* nearest neighbor */
	    k1 = floorf(x1); x1 -= k1;
	    k2 = floorf(x2); x2 -= k2;

	    switch(interp[0]) {
		case 'n': /* nearest neighbor */
		    if (x1 > 0.5) k1++;
		    if (x2 > 0.5) k2++;
		    if (k1 >=0 && k1 < n1 &&
			k2 >=0 && k2 < n2) {
			rotd[i2][i1] = orig[k2][k1];
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		case 'l': /* bilinear */
		    if (k1 >=0 && k1 < n1-1 &&
			k2 >=0 && k2 < n2-1) {
			rotd[i2][i1] = 
			    (1.-x1)*(1.-x2)*orig[k2][k1]   +
			    x1     *(1.-x2)*orig[k2][k1+1] +
			    (1.-x1)*x2     *orig[k2+1][k1] +
			    x1     *x2     *orig[k2+1][k1+1];
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		case 'c': /* cubic convolution */
		    /* !!! ADD CODE !!! */
		    break;
		default:
		    sf_error("Unknown interpolation %s",
			     interp);
		    break;
	    }
	}
    }

    /* write result */
    sf_floatwrite(rotd[0],n1*n2,out);

    exit(0);
}
