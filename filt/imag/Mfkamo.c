#include <math.h>
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nw, nx, ny, iw, ix, iy;
    float complex *oper;
    float dw,dx,dy, ow,ox,oy, w,x,y, x1,x2, h1,h2,f1,f2, maxe;
    float eps1,eps2,amp1,amp2,phase1,phase2,amp;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");

    if (!sf_histfloat (in,"o1",&ow)) sf_error("No o1= in input");
    if (!sf_histfloat (in,"d1",&dw)) sf_error("No d1= in input");
    if (!sf_histfloat (in,"o2",&ox)) sf_error("No o2= in input");
    if (!sf_histfloat (in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat (in,"o3",&oy)) sf_error("No o3= in input");
    if (!sf_histfloat (in,"d3",&dy)) sf_error("No d3= in input");

    if (!sf_getfloat("h1",&h1)) sf_error("Need h1=");
    if (!sf_getfloat("h2",&h2)) sf_error("Need h2=");
    if (!sf_getfloat("f1",&f1)) sf_error("Need f1=");
    if (!sf_getfloat("f2",&f2)) sf_error("Need f2=");

    if (!sf_getfloat("maxe",&maxe)) maxe=10.;

    f1 *= SF_PI/180.;
    f2 *= SF_PI/180.;

    oper = sf_complexalloc (nw);

    for (iy=0; iy < ny; iy++) {
	y = oy + iy*dy;
	for (ix=0; ix < nx; ix++) {
	    x = ox + ix*dx;
	    x1 = x*cos(f1) + y*sin(f1);
	    x2 = x*cos(f2) + y*sin(f2);
	    for (iw=0; iw < nw; iw++) {
		w = ow + iw*dw;
		if (fabsf (w) > FLT_EPSILON) {
		    eps1 = 2.*fabsf(x1*h1/w);
		    eps2 = 2.*fabsf(x2*h2/w);
		    if (eps1 <= maxe && eps2 <= maxe) {
			eps1 = sqrtf (1.+eps1*eps1);
			eps2 = sqrtf (1.+eps2*eps2);
                 
			amp1 = 1./eps1+eps1;
			amp2 = 1./eps2+eps2;
			phase1 = 1-eps1+logf(0.5*(1.+eps1));
			phase2 = 1-eps2+logf(0.5*(1.+eps2));

			amp = expf(0.5*(eps1-logf(amp1)+logf(amp2)-eps2));
			oper[iw] = amp*cexpf(I*SF_PI*(phase1-phase2)*w);
		    } else {
			oper[iw] = 0.;
		    }
		} else {
		    oper[iw] = 0.;
		}
	    }
	    sf_write (oper,sizeof(float complex),nw,out);
	}
    }

    exit (0);
}

