/* Post-stack 2-D v(z) time modeling/migration with curved-ray phase-shift.

Takes: < input.rsf > output.rsf
*/

#include <rsf.h>

#include "gazdagz.h"

int main (int argc, char *argv[])
{
    int nt;		/* number of time samples */
    int nz;		/* number of migrated time samples */
    int nx;		/* number of midpoints 	*/
    int ik,ix,it,iz;    /* loop counters 	*/
    int nxfft;	        /* fft size		*/
    int nk;		/* number of wave numbers */	
    
    float dt;		/* time sampling interval 	*/
    float dz;		/* migrated time sampling interval */
    float dk;	        /* wave number sampling interval */
    float k2;             /* wave number squared */
    float dx;		/* spatial sampling interval	*/
    float *vt, v0;	/* velocity v(t)		*/
    float **p,**q;	/* input, output data		*/

    float complex **cp,**cq; /* complex input,output	*/

    bool inv;             /* modeling or migration        */
    bool depth;           /* time or depth migration      */
    bool midpoint;        /* use midpoint rule            */
    float eps;            /* dip filter constant          */               

    sf_file in, out, vel;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* Stabilization parameter */

    if (!sf_histint(in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) 
	sf_error ("No d2= in input");

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");

	if (!sf_getint("nt",&nt)) sf_error ("nt= must be supplied");
	/* Length of time axis (for modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error ("dt= must be supplied");
        /* Sampling of time axis (for modeling) */

	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1= in input");
	if (NULL == sf_getstring("velocity")) {
	    if (!sf_getint("nz",&nz)) nz = nt;
	    /* Length of depth axis (for migration, if no velocity file) */
	    if (!sf_getfloat("dz",&dz)) dz = dt;
	    /* Sampling of depth axis (for migration, if no velocity file) */
	} else {
	    vel = sf_input("velocity");
	    if (!sf_histint(vel,"n1",&nz)) 
		sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) 
		sf_error ("No d1= in velocity");
	}
	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
    }

    vt = sf_floatalloc(nz);
    if (NULL == sf_getstring("velocity")) {
	/* file with velocity */
	if (!sf_getfloat("vel",&v0)) sf_error ("vel= must be supplied");
	/* Constant velocity (if no velocity file) */
	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	}
    } else {
	vel = sf_input("velocity");
	sf_floatread(vt,nz,vel);
	sf_fileclose(vel);
    }

    if (!sf_getbool("depth",&depth)) depth = false;
    /* if true, depth migration */

    if (!sf_getbool("midpoint",&midpoint)) midpoint = false;
    /* if true, use midpoint rule */

    /* vt -> 1/4 vt^2 */ 
    for (iz=0; iz < nz; iz++) {
	vt[iz] *= 0.25*vt[iz];
	if (depth) vt[iz] = 1./vt[iz];
    }

 
    /* determine wavenumber sampling (for real to complex FFT) */
    nxfft = sf_npfar(nx);
    nk = nxfft/2+1;
    dk = 2.0*SF_PI/(nxfft*dx);
	
    /* allocate space */
    p = sf_floatalloc2(nt,nxfft);
    q = sf_floatalloc2(nz,nxfft);
    cp = sf_complexalloc2(nt,nk);
    cq = sf_complexalloc2(nz,nk);

    if (inv) {
	sf_floatread(q[0],nz*nx,in);
    
	/* pad with zeros and Fourier transform x to k */
	for (ix=nx; ix<nxfft; ix++) {
	    for (iz=0; iz<nz; iz++) {
		q[ix][iz] = 0.0;
	    }
	}
    
	sf_pfa2rc(-1,2,nz,nxfft,q[0],cq[0]);
    } else {
	sf_floatread(p[0],nt*nx,in);
    
	/* pad with zeros and Fourier transform x to k */
	for (ix=nx; ix<nxfft; ix++) {
	    for (it=0; it<nt; it++) {
		p[ix][it] = 0.0;
	    }
	}
    
	sf_pfa2rc(-1,2,nt,nxfft,p[0],cp[0]);
    }

    gazdagz_init(eps,nt,dt,nz,dz,vt,depth,midpoint); /* ,gt); */

    /* migrate each wavenumber */
    for (ik=0; ik<nk; ik++) {
	k2 = ik*dk;
	k2 *= k2;
	gazdagz(inv,k2,cp[ik],cq[ik]);
    }	

    gazdagz_close();

    if (inv) {
	/* Fourier transform k to x (including FFT scaling) */
	sf_pfa2cr(1,2,nt,nxfft,cp[0],p[0]);
	for (ix=0; ix<nx; ix++) {
	    for (it=0; it<nt; it++) {
		p[ix][it] /= nxfft;
	    }
	}
	
	sf_floatwrite (p[0],nt*nx,out);
    } else {
	/* Fourier transform k to x (including FFT scaling) */
	sf_pfa2cr(1,2,nz,nxfft,cq[0],q[0]);
	for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
		q[ix][iz] /= nxfft;
	    }
	}

	sf_floatwrite (q[0],nz*nx,out);
    }


    exit (0);
}

/* 	$Id$	 */
