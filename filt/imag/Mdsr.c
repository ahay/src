#include <stdio.h>

#include <rsf.h>

#include "dsr.h"

int main (int argc, char **argv)
{
    int nt;               /* number of time samples */
    int nw;		/* number of frequencies */
    int nz;		/* number of migrated time samples */
    int nh;               /* number of offsets */
    int ik,im,jm,iz;      /* loop counters 	*/
    int nk,nm;		/* number of wave numbers */	

    float dt;             /* time sampling interval */
    float dw;		/* frequency sampling interval 	*/
    float w0;             /* first frequency */
    float t0;             /* first time */
    float dz;		/* migrated time sampling interval */
    float z0;             /* first migrated time */
    float dh;             /* offset sampling interval */
    float dk,dm;	        /* wave number sampling interval */
    float k,m;            /* wave number                  */
    float k0=0.,m0;       /* first offset wavenumber */
    float *vt, v0;	/* velocity v(t)		*/

    float complex **cp,*cq;	/* complex input,output		*/

    bool inv;              /* modeling or migration        */
    float eps;            /* dip filter constant          */               
    sf_file vel, in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    if (!sf_getfloat("eps",&eps)) eps = 0.01;

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");

	if (!sf_histint(in,"n2",&nk)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dk)) sf_error ("No d2= in input");

	if (!sf_getint("nt",&nt)) sf_error ("Need nt=");
	if (!sf_getfloat("dt",&dt)) sf_error ("Need dt=");
	if (!sf_getfloat("t0",&t0)) t0 = 0.;
	
	sf_putint(out,"nt",nt);
	sf_putfloat(out,"t0",t0);

	/* determine frequency sampling */
	nw = sf_npfa(nt);
	dw = 2.0*SF_PI/(nw*dt);
	w0 = -SF_PI/dt;

	sf_putint(out,"n1",nw);
	sf_putfloat(out,"d1",dw);
	sf_putfloat(out,"o1",w0);
	
	if (!sf_getint("nh",&nh)) sf_error ("Need nh=");
	if (!sf_getfloat("dh",&dh)) sf_error ("Need dh=");
    
	/* determine wavenumber sampling, pad by 2 */
	nm = nh*2;
	nm = sf_npfao(nm,nm*2);
	dm = 2.0*SF_PI/(nm*dh);
	m0 = -SF_PI/dh;

	sf_putint(out,"n2",nm);
	sf_putfloat(out,"d2",dm);
	sf_putfloat(out,"o2",m0);

	/* ny? */
	sf_putint(out,"ny",nh);
	sf_putfloat(out,"dy",dh);

	sf_putint(out,"n3",nk);
	sf_putfloat(out,"d3",dk);
	sf_putfloat(out,"o3",k0);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nw)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dw)) sf_error ("No d1= in input");
	if (!sf_histfloat(in,"o1",&w0)) sf_error ("No o1= in input");

	if (!sf_histint(in,"n2",&nm)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dm)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&m0)) sf_error ("No o2= in input");
	
	if (!sf_histint(in,"n3",&nk)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dk)) sf_error ("No d3= in input");

	if (NULL == sf_getstring("velocity")) {
	    vel = NULL;
	    if (!sf_getint("nz",&nz)) sf_error ("Need nz=");
	    if (!sf_getfloat("dz",&dz)) sf_error ("Need dz=");
	    if (!sf_getfloat("z0",&z0)) z0 = 0.;
	} else {
	    vel = sf_input ("velocity");
	    if (!sf_histint(vel,"n1",&nz)) sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) sf_error ("No d1= in velocity");
	    if (!sf_histfloat(vel,"o1",&z0)) z0 = 0.; 
	}
	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",z0);

	sf_putint(out,"n2",nk);
	sf_putfloat(out,"d2",dk);
	sf_putfloat(out,"o2",k0);

	sf_putint(out,"n3",1);
    }

    vt     = sf_floatalloc(nz);

    if (NULL == vel) {
	if (!sf_getfloat("vel",&v0)) sf_error ("Need vel=");
	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	}
    } else {
	sf_read(vt,sizeof(float),nz,vel);
    }

    /* allocate space */
    cp = sf_complexalloc2(nw,nm);
    cq = sf_complexalloc(nz);

    /* migrate each wavenumber */
    for (ik=0; ik<nk; ik++) {
	sf_warning("wavenumber %d of %d",ik+1,nk);
	
	k = ik*dk;

	if (inv) { /* modeling */
	    sf_read(cq,sizeof(float complex),nz,in);
	} else {
	    for (iz=0; iz<nz; iz++) {
		cq[iz] = 0.0;
	    }
	    sf_read(cp[0],sizeof(float complex),nw*nm,in);
	}

	for (im=0; im<nm; im++) {
	    jm = (im < nm/2)? im + nm/2: im - nm/2;
	    m = m0 + jm*dm;
      
	    dsr(inv,eps,k,m,nw,dw,w0,nz,dz,vt,cp[im],cq);
	}
 
	if (inv) {
	    sf_write(cp[0],sizeof(float complex),nw*nm,out);
	} else {
	    sf_write(cq,sizeof(float complex),nz,out);
	}
    }

    exit (0);
}

