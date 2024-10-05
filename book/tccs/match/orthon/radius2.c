/* smoothing radii */

#include <rsf.h>

int main (int argc, char* argv[])
{
    int n1, n1f, n2, n2f, i, n12, n12f;
    float *rect, *fr, maxrad, c, *rad, *rad_low, *rad_high;
    sf_file in, out, freq, low, high;
	
    sf_init (argc,argv);

    in = sf_input("in");
    freq = sf_input("freq");
    out = sf_output("out");
    low = sf_output("low");
    high = sf_output("high");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(freq,"n1",&n1f)) sf_error("No n1= in frequency difference.");

    n2 = sf_leftsize(in,1);    
    n2f = sf_leftsize(freq,1);    
    
    n12 = n1*n2;
    n12f = n1f*n2f;

    if (n1 != n1f) sf_error("Need matching n1");
    if (n2 != n2f) sf_error("Need matching n2");

    if (!sf_getfloat("c",&c)) c=1.;
    if (!sf_getfloat("maxrad",&maxrad)) maxrad=1000.;

    rect = sf_floatalloc(n12);
    sf_floatread(rect,n12,in);

    fr = sf_floatalloc(n12f);
    sf_floatread(fr,n12,freq);

    rad = sf_floatalloc(n12);
    rad_low = sf_floatalloc(n12);
    rad_high = sf_floatalloc(n12);

    /* constraint conditions: [-maxrad, -1] U [1, maxrad] */
    for (i=0; i < n12; i++) {

        /* update radius */
	    rad[i] = rect[i]+c*fr[i];

        /* set maximum allowed radius */
        if (rad[i] > maxrad) rad[i] = maxrad;

        /* low radius */
        if (rad[i] < 0){
           rad_high[i] = 1;
           rad_low[i] = -1 * rad[i];
           if (rad_low[i] < 1) rad_low[i] = 1;

        /* high radius */
        }else{
           rad_low[i] = 1;
           rad_high[i] = rad[i];
           if (rad_high[i] < 1) rad_high[i] = 1;
        }
    }

    sf_floatwrite(rad,n12,out);
    sf_floatwrite(rad_low,n12,low);
    sf_floatwrite(rad_high,n12,high);
    exit(0);
}