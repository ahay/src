/* Velocity continuation with semblance computation.

Takes < input.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

#include "stretch.h"
#include "cosft.h"

int main(int argc, char* argv[])
{
    map str, istr;
    int i1,i2, n1,n2,n3, ix,iv,ih, ib, ie, nb, nx,nv,nh, nw;
    float d1,o1,d2,o2, eps, w,x,k, v0,v2,v,v1,dv, dx, h0,dh,h, num, den, t;
    float *trace, *strace, ***stack, ***stack2, ***cont, **image;
    float complex *ctrace, *ctrace0;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input");

    if (!sf_getint("nb",&nb)) nb=2;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getint("pad",&n2)) n2=n1;
    if (!sf_getint("pad2",&n3)) n3=n2;

    n3 = sf_npfar(n3);
    nw = n3/2+1;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;  
    o2 = o1*o1;
 
    if(!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    d2 = o1+(n1-1)*d1;
    d2 = (d2*d2 - o2)/(n2-1);

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    if (!sf_getfloat("v0",&v0) && 
	!sf_histfloat(in,"v0",&v0)) sf_error("Need v0=");

    if(!sf_histfloat(in,"o3",&h0)) sf_error("No o2= in input");
    if(!sf_histfloat(in,"d3",&dh)) sf_error("No d2= in input");
    if(!sf_histfloat(in,"d2",&dx)) sf_error("No d3= in input");
    
    sf_putfloat(out,"o3",v0+dv);
    sf_putfloat(out,"d3",dv);
    sf_putint(out,"n3",nv);

    sf_putstring(out,"label3","Velocity (km/s)");

    dx = 2.*SF_PI/(sf_npfar(2*(nx-1))*dx);

    stack = sf_floatalloc3(n1,nx,nv);
    stack2 = sf_floatalloc3(n1,nx,nv);
    cont = sf_floatalloc3(n1,nx,nv);
    image = sf_floatalloc2(n1,nx);
    trace = sf_floatalloc(n1);
    strace = sf_floatalloc(n3);
    ctrace = sf_complexalloc(nw);
    ctrace0 = sf_complexalloc(nw);

    for (i1=0; i1 < n1; i1++) {
	t = o1+i1*d1;
	trace[i1] = t*t;
    }

    str = stretch_init (n2, o2, d2, n1, eps);
    stretch_define (str, trace);

    for (i2=0; i2 < n2; i2++) {
	t = o2+i2*d2;
	strace[i2] = sqrtf(t);
    }

    istr = stretch_init (n1, o1, d1, n2, eps);
    stretch_define (istr, strace);

    for (i1=0; i1 < n1*nx*nv; i1++) {
	stack[0][0][i1] = 0.;
	stack2[0][0][i1] = 0.;
    }

    cosft_init(nx);

    for (ih=0; ih < nh; ih++) {
	sf_warning("offset %d of %d",ih+1,nh);

	h = h0 + ih*dh;
	h *= h;

	sf_read(image[0],sizeof(float),n1*nx,in);

	for (i1=0; i1 < n1; i1++) {
	    cosft_frw(image[0],i1,n1);
	}

	for (ix=0; ix < nx; ix++) {
	    x = ix*dx; 
	    x *= x;

	    k = x * 0.25 * 0.25 * 0.5;

	    stretch_apply (str, image[ix], strace);
	    
	    for (i2=n2; i2 < n3; i2++) {
		strace[i2] = 0.;
	    }
		
	    sf_pfarc(1,n3,strace,ctrace0);

	    for (iv=0; iv < nv; iv++) {
		v = v0 + (iv+1)* dv;

		v1 = h * (1./(v*v) - 1./(v0*v0)) * 8.;
		v2 = k * ((v0*v0) - (v*v));

		ctrace[0]=0.; /* dc */

		for (i2=1; i2 < nw; i2++) {
		    w = i2*SF_PI/(d2*n3);
 
		    ctrace[i2] = ctrace0[i2] * cexpf(-I*(v2/w+(v1-o2)*w));
		} /* w */

		sf_pfacr(-1,n3,ctrace,strace);

		stretch_apply (istr, strace, cont[iv][ix]);    
	    } /* v */
	} /* x */

	for (iv=0; iv < nv; iv++) {
	    for (i1=0; i1 < n1; i1++) {
		cosft_inv(cont[0][0],i1+iv*nx*n1,n1);
	    }
	}

	for (iv=0; iv < nv; iv++) {
	    for (ix=0; ix < nx; ix++) {
		for (i1=0; i1 < n1; i1++) {	
		    t = cont[iv][ix][i1];
		    stack[iv][ix][i1] += t;
		    stack2[iv][ix][i1] += t*t;
		} /* i1 */
	    } /* x */
        } /* v */
    } /* h */

    for (iv=0; iv < nv; iv++) {
	for (ix=0; ix < nx; ix++) {
	    for (i1=0; i1 < n1; i1++) {
		ib = i1-nb > 0? i1-nb: 0;
		ie = i1+nb+1 < n1? i1+nb+1: n1;

		num = 0.;
		den = 0.;

		for (i2=ib; i2 < ie; i2++) {
		    t = stack[iv][ix][i2];
		    num += t*t;
		    den += stack2[iv][ix][i2];
		}

		trace[i1] = den > 0.? num/den: 0.;
	    }
	    sf_write (trace, sizeof(float), n1, out);
	}
    }

    exit(0);
}


