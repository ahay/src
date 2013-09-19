/* Analytical first-arrival traveltimes. */
#include <math.h>
#include <rsf.h>

int main (int argc, char* argv[])
{
    char *type;
    int n1, n2, i1, i2;
    float d1,d2, g1,g2, s,v0, x1,x2,gsq,g,s2,z,d;
    float *time;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* Get grid dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1=");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2=");
    
    if (!sf_getfloat("g1",&g1)) g1 = 0.;    
    /* vertical gradient */
    if (!sf_getfloat("g2",&g2)) g2 = 0.;    
    /* horizontal gradient */
    gsq = g1*g1+g2*g2;
    g = sqrtf(gsq);

    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* initial velocity or slowness squared */

    if (!sf_getfloat("s",&s)) s = 0.0;
    /* shot location at the surface */
 
    if (NULL == (type = sf_getstring("case"))) type="c";
    /* case of velocity distribution */

    if (0.0 == g1 && 0.0 == g2) type="const";
    
    time = sf_floatalloc(n1);

    for (i2 = 0; i2 < n2; i2++) {
	x2 = i2*d2;
	for (i1 = 0; i1 < n1; i1++) {
	    x1 = i1*d1;
	    d = x1*x1+(x2-s)*(x2-s);
	    switch (type[0]) {
		case 's':
		    /* slowness squared */
		    s2 = v0+g1*x1+g2*x2;
		    z = 2.0*d/(s2+sqrtf(s2*s2-gsq*d));
		    time[i1] = (s2-gsq*z/6.0)*sqrtf(z);
		    break;
		case 'v': 
		    /* velocity */
		    s2 = 2.0*v0*(v0+g1*x1+g2*x2);
		    /* !!! CHANGE BELOW !!! */
		    time[i1] = hypotf(x2-s,x1)/v0;
		    break;
		case 'c': /* constant velocity */
		default:
		    time[i1] = hypotf(x2-s,x1)/v0;
		    break;
	    }
	}
	sf_floatwrite(time,n1,out);
    }

    exit (0);
}
