/* DMO and stack by finite-difference offset continuation.

Takes: < nmod.rsf > stack.rsf
*/
#include <math.h>

#include <rsf.h>

#include "ctridiagonal.h"

int main(int argc, char* argv[])
{
    ctris slv;
    bool stack;
    int nw,nh,nx,iw,ix,ih,n1,n2,ns;
    float w0,dw, h0,dh,dx,w,w2,h,h2;
    float complex diag,diag2, c1,c2,offd,offd2;
    float complex *in, *out, **dat;
    sf_file cmp, stk;

    sf_init (argc,argv);
    cmp = sf_input("in");
    stk = sf_output("out");

    if (SF_COMPLEX != sf_gettype(cmp)) sf_error("Need complex input");
    if (!sf_histint(cmp,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(cmp,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
    dh /= dx;
    h0 /= dx;
    if (!sf_histfloat(cmp,"d3",&dw)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&w0)) sf_error("No o3= in input");

    if (!sf_getbool("stack",&stack)) stack=true;
    /* if n, don't stack, output everything */

    in = sf_complexalloc(nx);
    out = sf_complexalloc(nx);
    dat = sf_complexalloc2(nx,nh);

    slv = ctridiagonal_init (nx);

    if (dh > 0.) {
	n1 = nh-1; 
	ns = -1;
	n2 = -0.5 - h0/dh;		    
    } else {
	n1 = 0; 
	ns = 1;
	n2 = 0.5 - h0/dh;
    }

    if (stack) {
	sf_putint(stk,"n3",1);
	sf_putint(stk,"n2",nw);
	sf_putfloat(stk,"o2",w0);
	sf_putfloat(stk,"d2",dw);
    } else {
	sf_putint(stk,"n2",(n2-n1+ns)/ns);
	sf_putfloat(stk,"o2",h0+(n2-n1)*dh);
    }

    for (iw=0; iw < nw; iw++) {
	sf_warning("frequency %d of %d",iw+1,nw);

	w = -2.*SF_PI*(w0 + iw*dw); 
	w2 = w*w;

	sf_read(dat[0],sizeof(float complex),nx*nh,cmp);

	for (ix=0; ix < nx; ix++) {
	    out[ix] = 0.;
	}

	if (fabsf(w) < dw) {	    
	    if (stack) {
		sf_write(out,sizeof(float complex),nx,stk);
	    } else {
		for (ih=n1; ih != n2; ih += ns) {
		    sf_write(out,sizeof(float complex),nx,stk);
		}
	    }
	    continue;
	}

	c1 = 3.*(9. + w2 + 4.*w*I)/(w2*(3. - w*I));
	c2 = 3.*(w2 - 27. + 8.*w*I)/(w2*(3. - w*I));

	for (ih=n1; ih != n2; ih += ns) {
	    for (ix=0; ix < nx; ix++) {
		if (ih >=0 && ih < nh) {		    
		    in[ix] = stack? dat[ih][ix] + out[ix]: dat[ih][ix];
		} else {
		    in[ix] = out[ix];
		}
	    }

	    h = h0 + ih*dh; 
	    h2 = h+dh*ns; 
	    h = h*h; 
	    h2 = h2*h2;

	    offd  = 1. - c1*h2 + c2*h;
	    offd2 = 1. - c1*h  + c2*h2;
	    diag  = 12. - 2.*offd;
	    diag2 = 12. - 2.*offd2;

	    ctridiagonal_const_define (slv,diag2,offd2);

	    out[0] = diag * in[0] + offd * in[1];
	    for (ix=1; ix < nx-1; ix++) {
		out[ix] = diag * in[ix] + offd * (in[ix+1] +in[ix-1]);
	    }
	    out[nx-1] = diag * in[nx-1] + offd * in[nx-2];

	    ctridiagonal_solve (slv,out);

	    if (!stack) sf_write (out,sizeof(float complex),nx,stk);
	}
	if (stack) sf_write (out,sizeof(float complex),nx,stk);
    }

    sf_close();
    exit(0);
}



