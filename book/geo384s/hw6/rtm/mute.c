/* Mute a triangle region */
#include <rsf.h>

int main(int argc, char* argv[]) {
    int n1, n2, n3, i1, i2, i3;
    float t0,v0, dt,dv, t1, v1, slope,tmax,band, t,taper, *data;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&v0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&dv)) sf_error("No d2= in input");

    if (!sf_getfloat("t1",&t1)) t1=t0; // start time
    if (!sf_getfloat("v1",&v1)) v1=v0+(n2-1)*dv; // end velocity
    if (!sf_getfloat("band",&band)) band=20*dt; // start time
    slope = (t0+(n1-1)*dt-t1)/(v1-v0);

    data = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
        for (i2=0; i2 < n2; i2++) {
            tmax = t1 + i2*dv*slope;
            sf_floatread (data,n1,in);

            for (i1=0; i1 < n1; i1++) {
                t = t0+i1*dt;

                if (t > tmax) {
                    data[i1] = 0.0f;
                } else if (t > tmax-band) {
                    taper = sinf(0.5 * SF_PI * (t-tmax)/band);
                    data[i1] *= taper*taper;
                }
            }
            sf_floatwrite (data,n1,out);
        }
    }

    exit(0);
}


