/* 2-D post-stack modeling/migration with split step.

Takes: < input.rsf > output.rsf
*/

#include <rsf.h>

#include "split1.h"

int main (int argc, char *argv[])
{
    int nt;		/* number of time samples */
    int nz;		/* number of migrated time samples */
    int nx;		/* number of midpoints 	*/
    int ix,it,iz;         /* loop counters 	*/
    int ntfft;	        /* fft size		*/
    int nw;		/* number of frequencies */	

    float dt;		/* time sampling interval 	*/
    float dz;		/* migrated time sampling interval */
    float dw;	        /* frequency sampling interval */
    float dx;		/* spatial sampling interval	*/
    float **vt, *v, v0;	/* velocities		*/
    float **p,**q;	/* input, output data		*/

    float complex **cp;	   /* complex input		*/

    bool inv;             /* modeling or migration        */
    bool depth;           /* depth or time                */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant          */  
    sf_file in, out, vel;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool("depth",&depth)) depth = false;
    /* depth or time migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stability parameter */

    if (!sf_histint(in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");

    if (NULL == sf_getstring("velocity")) {
	/* velocity file */
	if (!sf_getfloat("vel",&v0)) sf_error ("Need vel=");
	/* constant velocity (if no velocity file) */
	vel = NULL;	
    } else {
	vel = sf_input("velocity");
    }

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");
	if (!sf_getint("nt",&nt)) sf_error ("Need nt=");
	/* Length of time axis (for modeling) */ 
	if (!sf_getfloat("dt",&dt)) sf_error ("Need dt=");
	/* Time sampling (for modeling) */
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1= in input");
	if (NULL == vel) {
	    if (!sf_getint("nz",&nz)) {
		/* number of steps in depth 
		   (for constant-velocity depth migration) */
		if (depth) sf_error ("Need nz=");
		nz = nt;
	    }
	    if (!sf_getfloat("dz",&dz)) {
		/* sampling in depth 
		   (for constant-velocity depth migration) */
		if (depth) sf_error ("Need dz=");
		dz = dt*v0;
	    }
	} else {
	    if (!sf_histint(vel,"n1",&nz)) 
		sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) 
		sf_error ("No d1= in velocity");
	}
	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
    }
    
    vt = sf_floatalloc2(nz,nx);
    v = sf_floatalloc(nz);

    if (NULL == vel) {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		vt[ix][iz] = v0;
	    }
	}
    } else {
	sf_read(vt[0],sizeof(float),nx*nz,vel);
	sf_fileclose(vel);
    }
    for (iz=0; iz < nz; iz++) {
	v[iz] = 0.;
	for (ix=0; ix < nx; ix++) {
	    vt[ix][iz] = 4./(vt[ix][iz]*vt[ix][iz]);
	    v[iz] += vt[ix][iz];
	}
	v[iz] /= nx;
    }

    /* determine wavenumber sampling (for real to complex FFT) */
    ntfft = sf_npfaro(nt,nt*2);
    nw = ntfft/2+1;
    dw = 2.0*SF_PI/(ntfft*dt);
	
    /* allocate space */
    p = sf_floatalloc2(ntfft,nx);
    q = sf_floatalloc2(nz,nx);
    cp = sf_complexalloc2(nw,nx);

    for (ix=0; ix<nx; ix++) {
	if (inv) {
	    sf_read(q[ix],sizeof(float),nz,in);
	} else {    
	    sf_read(p[ix],sizeof(float),nt,in);

	    /* pad with zeros and Fourier transform t to w */
	    for (it=nt; it<ntfft; it++) {
		p[ix][it] = 0.0;
	    }
	    
	    sf_pfarc(-1,ntfft,p[ix],cp[ix]);
	    for (it=0; it<nw; it++)
		cp[ix][it] /= ntfft;
	}
    }

    split1 (verb, inv, eps,  
	    nw, dw, 
	    nz, dz, 
	    nx, dx,
	    vt, v,
	    cp, q);

    for (ix=0; ix<nx; ix++) {
	if (inv) {
	    /* Fourier transform w to t (including FFT scaling) */
	    sf_pfacr(1,ntfft,cp[ix],p[ix]);
	    for (it=0; it<nt; it++)
		p[ix][it] /= ntfft;
	    sf_write (p[ix],sizeof(float),nt,out);
	} else {
	    sf_write (q[ix],sizeof(float),nz,out);
	}
    }
    
    sf_close();
    exit (0);
}

/* 	$Id: Msstep1.c,v 1.4 2004/03/22 05:43:24 fomels Exp $	 */
