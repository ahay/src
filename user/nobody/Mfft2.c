/* 2-D Fast Fourier Transform.
*/
/*
Copyright (C) 2004 University of Texas at Austin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

int main (int argc, char **argv)
{
    int nt, nx, n3;	/* dimensions */
    int ik,ix,it,i3;      /* loop counters 	*/
    int nxfft;	        /* fft size		*/
    int nk;		/* number of wavenumbers */	
    int nw;               /* number of frequencies */

    float dx;		/* space sampling interval */
    float dk;	        /* wavenumber sampling interval */
    float dt;		/* time sampling interval	*/
    float dw;             /* frequency sampling interval */
    float w0;             /* starting frequency */
    float t0;             /* starting time */
    float x0;             /* staring space */
    float k0=0.;          /* starting wavenumber */

    float **p;	        /* time-space */
    float complex **cp;	        /* frequency-wavenumber */
    float complex *cq;           /* frequency */

    bool inv;              /* forward or inverse */
    bool both;             /* both coordinates or second only */
  
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* if y, perform inverse transform */
    if (!sf_getbool("both",&both)) both = false;
    /* if y, transform both axes */

    n3 = sf_leftsize(in,2);

    if (inv) { 
	if (!sf_histint(in,"n1",&nw)) sf_error ("No n1 in input");
	if (!sf_histfloat(in,"d1",&dw)) sf_error ("No d1 in input");
	if (!sf_histint(in,"n2",&nk)) sf_error ("No n2 in input");
	if (!sf_histfloat(in,"d2",&dk)) sf_error ("No d2 in input");

	if (!sf_histint(in,"nx",&nx)) sf_error ("No nx in input");
	if (!sf_histfloat(in,"x0",&x0)) x0 = 0.; 

	if (both) {
	    if (!sf_histint(in,"nt",&nt)) sf_error ("No nt in input");
	    if (!sf_histfloat(in,"t0",&t0)) t0 = 0.; 

	    dt = 1./(nw*dw);
      
	    sf_putint(out,"n1",nt);
	    sf_putfloat(out,"d1",dt);
	    sf_putfloat(out,"o1",t0);
	} else {
	    nt = nw;
	}

	nxfft = 2*(nk-1);
	dx = 1./(nxfft*dk);

	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);
	
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype(out,SF_FLOAT);
    } else { 
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1 in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1 in input");
	if (!sf_histfloat(in,"o1",&t0)) t0 = 0.; 
	if (!sf_histint(in,"n2",&nx)) sf_error ("No n2 in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2 in input");
	if (!sf_histfloat(in,"o2",&x0)) x0 = 0.;

	if (both) {
	    sf_putint (out,"nt",nt);
	    sf_putfloat (out,"t0",t0);

	    /* determine frequency sampling */
	    nw = sf_npfa(nt);
	    dw = 1./(nw*dt);
	    w0 = -0.5/dt;

	    sf_putint(out,"n1",nw);
	    sf_putfloat(out,"d1",dw);
	    sf_putfloat(out,"o1",w0);
	} else {
	    nw = nt;
	}

	sf_putint(out,"nx",nx);
	sf_putfloat(out,"x0",x0);

	/* determine wavenumber sampling (for real to complex FFT) */
	nxfft = nx*2;
	nxfft = sf_npfaro(nxfft,2*nxfft);
	nk = nxfft/2+1;
	dk = 1./(nxfft*dx);

	sf_putint(out,"n2",nk);
	sf_putfloat(out,"d2",dk);
	sf_putfloat(out,"o2",k0);

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
	sf_settype(out,SF_COMPLEX);
    }

    p = sf_floatalloc2(nt,nxfft);
    cp = sf_complexalloc2(nt,nk);
    cq = sf_complexalloc(nw);

    for (i3=0; i3<n3; i3++) {
	if (inv) {
	    for (ik=0; ik<nk; ik++) {
		sf_complexread(cq,nw,in);
		
		if (ik==0) {
		    for (it=0; it<nt; it++) {
			cp[ik][it] = 0.;
		    }
		} else if (both) {
		    /* Fourier transform w to t */
		    sf_pfacc(-1,nw,cq);
		    for (it=0; it<nt; it++) {
			cp[ik][it] = (it%2 ? -cq[it] : cq[it]);
		    }
		} else {
		    for (it=0; it<nt; it++) {
			cp[ik][it] = cq[it];
		    }
		}
	    }

	    /* Fourier transform k to x */
	    sf_pfa2cr(1,2,nt,nxfft,cp[0],p[0]);

	    /* FFT scaling */
	    for (ix=0; ix<nx; ix++) {
		for (it=0; it<nt; it++) {
		    p[ix][it] /= nxfft;
		}
	    }

	    sf_floatwrite(p[0],nt*nx,out);
	} else { /* forward */
	    sf_floatread(p[0],nt*nx,in);
      
	    /* pad with zeros */
	    for (ix=nx; ix<nxfft; ix++) {
		for (it=0; it<nt; it++) {
		    p[ix][it] = 0.0;
		}
	    }
    
	    /* Fourier transform x to k */
	    sf_pfa2rc(-1,2,nt,nxfft,p[0],cp[0]);

	    for(ik=0; ik<nk; ik++) {
		if (ik==0) {
		    for (it=0; it<nw; it++) {
			cq[it] = 0.;
		    }
		} else if (both) {
		    /* Fourier transform t to w, with w centered */
		    for (it=0; it<nt; it++) {
			/* include FFT scaling */
			cq[it] = (it%2 ? -cp[ik][it] : cp[ik][it])/nw;
		    }
		    /* Pad with zeros */
		    for (it=nt; it<nw; it++) {
			cq[it] = 0.;
		    }  
		    sf_pfacc(1,nw,cq);
		} else {
		    for (it=0; it<nt; it++) {
			cq[it] = cp[ik][it];
		    }
		    /* Pad with zeros */
		    for (it=nt; it<nw; it++) {
			cq[it] = 0.;
		    }  
		}
	
		sf_complexwrite(cq,nw,out); 
	    }
	}
    }

    exit (0);
}

/* 	$Id$	 */
