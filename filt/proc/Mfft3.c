#include <rsf.h>

int main (int argc, char **argv)
{
    int n1, nx, n3;	/* dimensions */
    int i1, ix, i3;       /* loop counters 	*/
    int nk;		/* number of wavenumbers */	

    float dx;		/* space sampling interval */
    float dk;	        /* wavenumber sampling interval */
    float x0;             /* staring space */
    float k0;             /* starting wavenumber */

    float complex **cp;	        /* frequency-wavenumber */

    bool inv;              /* forward or inverse */
    
    sf_file in, out;

    /* hook up getpar to handle the parameters */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;

    if (!sf_histint(in,"n1",&n1)) n1=1; 
    n3 = sf_leftsize(in,2);
    
    if (SF_COMPLEX != sf_gettype(in)) sf_error ("Need complex input");

    if (inv) { 
	if (!sf_histint(in,"n2",&nk)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dk)) sf_error ("No d2= in input");

	if (!sf_histint(in,"ny",&nx)) sf_error ("No nx= in input");
	if (!sf_histfloat(in,"y0",&x0)) x0 = 0.; 

	dx = 2.0*SF_PI/(nk*dk);

	sf_putint (out,"n2",nx);
	sf_putfloat (out,"d2",dx);
	sf_putfloat (out,"o2",x0);
    } else { 
	if (!sf_histint(in,"n2",&nx)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&x0)) x0 = 0.;

	sf_putint(out,"ny",nx);
	sf_putfloat(out,"y0",x0);

	/* determine wavenumber sampling, pad by 2 */
	nk = nx*2;
	nk = sf_npfao(nk,nk*2);
	dk = 2.0*SF_PI/(nk*dx);
	k0 = -SF_PI/dx;

	sf_putint (out,"n2",nk);
	sf_putfloat (out,"d2",dk);
	sf_putfloat (out,"o2",k0);
    }

    cp = sf_complexalloc2(n1,nk);

    for (i3=0; i3<n3; i3++) {
	if (inv) {
	    sf_read(cp[0],sizeof(float complex),n1*nk,in);
      
	    /* Fourier transform k to x */
	    sf_pfa2cc(1,2,n1,nk,cp[0]);

	    /* FFT scaling */
	    for (ix=0; ix<nx; ix++) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1] /= nk;
		}
	    }
      
	    sf_write(cp[0],sizeof(float complex),n1*nx,out);
	} else {
	    sf_read(cp[0],sizeof(float complex),n1*nx,in);
      
	    /* pad with zeros */
	    for (ix=nx; ix<nk; ix++) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1] = 0.;
		}
	    }
    
	    /* Fourier transform x to k */
	    sf_pfa2cc(-1,2,n1,nk,cp[0]);
      
	    ix = nk/2+1; /* oddball negative nyquist */
	    for (i1=0; i1<n1; i1++) {
		cp[ix][i1] = 0.;
	    }

	    sf_write(cp[0],sizeof(float complex),n1*nk,out);
	}
    }

    exit (0);
}

