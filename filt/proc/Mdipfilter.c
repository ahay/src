/*
#	pass=1	Pass dips in band
#	    =0	Reject dips in band
#
#	angle=0 Define band in terms of velocities
#	     =1 Define band in terms of angles
#
#	dim=2	2D dip-filter
#          =3	3D dip-filter
#
#	taper=1 linear taper
#            =2 cosine taper
#
#	Velocity/Dip band defined by:
#
#	    v2_____v3
#	     /     \
#	    /       \
#	 __/         \___
#	  v1         v4
#
# AUTHOR: James Rickett - Oct 97
*/
#include <float.h>
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1,n2,n3, nw, n2pad, nk, i3, i2, i1;
    float d1,d2, dk,k, w0,dw,w, vel,v;
    float ang1, ang2, ang3, ang4, v1, v2, v3, v4, factor, **data;
    float complex **cdata, *ctrace;
    bool angle, pass;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

/* determine frequency sampling (for complex FFT) */
    nw = sf_npfa(n1);
    dw = 2.0*SF_PI/(nw*d1);
    w0 = -SF_PI/d1;

/* determine wavenumber sampling (for real to complex FFT) */
    n2pad = sf_npfar(n2);
    nk = n2pad/2+1;
    dk = 2.0*SF_PI/(n2pad*d2);

    data = sf_floatalloc2(n1,n2pad);
    cdata = sf_complexalloc2(n1,nk);
    ctrace = sf_complexalloc(nw);

    if (!sf_getbool("angle",&angle)) angle=false;

    if (angle) {	
	if (!sf_getfloat("v",&v)) v=-1.;
	if (v < 0.) v=(d2*(n2-1))/(d1*(n1-1));

	if (!sf_getfloat("ang1",&ang1)) ang1=-50.;
	if (!sf_getfloat("ang2",&ang2)) ang2=-45.;
	if (!sf_getfloat("ang3",&ang3)) ang3=45.;
	if (!sf_getfloat("ang4",&ang4)) ang4=50.;
	v1=tanf(SF_PI*ang1/180.)*v;
	v2=tanf(SF_PI*ang2/180.)*v;
	v3=tanf(SF_PI*ang3/180.)*v;
	v4=tanf(SF_PI*ang4/180.)*v;
    } else {
	if (!sf_getfloat("v1",&v1)) v1=0.;
	if (!sf_getfloat("v2",&v2)) v2=0.1;
	if (!sf_getfloat("v3",&v3)) v3=99999.;
	if (!sf_getfloat("v4",&v4)) v4=999999.;
    }
    
    if ((v1>=v2)||(v2>=v3)||(v3>=v4)) 
	sf_error("Need v1 < v2 < v3 < v4, got %g ? %g ? %g ? %g",v1,v2,v3,v4); 

    if (!sf_getbool("pass",&pass)) pass=true;

    for (i3=0; i3 < n3; i3++) { 
	sf_read(data[0],sizeof(float),n1*n2,in);

	/* pad */ 
	for (i2=n2; i2 < n2pad; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		data[i2][i1]=0.;
	    }
	}

	/* Fourier transform x to k */
	sf_pfa2rc(-1,2,n1,n2pad,data[0],cdata[0]);

	for (i2=0; i2 < nk; i2++) {
	    /* Fourier transform t to w, with w centered */
	    for (i1=0; i1 < n1; i1++) {
		/* include FFT scaling */
		ctrace[i1] = (i1%2 ? -cdata[i2][i1] : cdata[i2][i1])/nw;
	    }	    
	    for (i1=n1; i1 < nw; i1++) {
		ctrace[i1] = 0.;
	    }
	    sf_pfacc(1,nw,ctrace);
	    
	    for (i1=0; i1 < nw; i1++) {
		w = w0+i1*dw;		    
		vel = w/(k+FLT_EPSILON); 
		
		if ((vel>=v1) && (vel<=v2)) {
		    factor=1.-sinf(0.5*SF_PI*(vel-v1)/(v2-v1));
		} else if ((vel>=v2) && (vel<=v3)) {
		    factor=0.;
		} else if ((vel>=v3) && (vel<=v4)) {
		    factor=sinf(0.5*SF_PI*(vel-v3)/(v4-v3));
		} else {
		    factor=1.;
		}
		
		if (pass) factor=1.-factor;
		
		ctrace[i1] *= factor;
	    }

	    /* Fourier transform w to t */
	    sf_pfacc(-1,nw,ctrace);

	    for (i1=0; i1 < n1; i1++) {
		cdata[i2][i1] = (i1%2 ? -ctrace[i1] : ctrace[i1]);
	    }
	}

	/* Fourier transform k to x */
	sf_pfa2cr(1,2,n1,n2pad,cdata[0],data[0]);

	/* FFT scaling */
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		data[i2][i1] /= n2pad;
	    }
	}

	sf_write(data[0],sizeof(float),n1*n2,out);
    } /* loop over n3 */

    exit(0);
}
