/* Prestack 2-D v(z) time modeling/migration by DSR.
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
    /* If y, modeling; If n, migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* Stabilization parameter */

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");

	if (!sf_histint(in,"n2",&nk)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dk)) sf_error ("No d2= in input");

	if (!sf_getint("nt",&nt)) sf_error ("Need nt=");
	/* Length of time axis (for modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error ("Need dt=");
	/* Sampling of time axis (for modeling) */
	if (!sf_getfloat("t0",&t0)) t0 = 0.;
	/* Origin of time axis (for modeling) */
	
	sf_putint(out,"nt",nt);
	sf_putfloat(out,"t0",t0);

	/* determine frequency sampling */
	nw = nt;
	dw = 1./(nw*dt);
	w0 = -0.5/dt;

	sf_putint(out,"n1",nw);
	sf_putfloat(out,"d1",dw);
	sf_putfloat(out,"o1",w0);

	if (!sf_getint("nh",&nh)) sf_error ("Need nh=");
	/* Number of offsets (for modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error ("Need dh=");
	/* Offset sampling (for modeling) */

	/* determine wavenumber sampling, pad by 2 */
	nm = nh*2;
	nm *= 2;
	dm = 1./(nm*dh);
	m0 = -0.5/dh;

	sf_putint(out,"n2",nm);
	sf_putfloat(out,"d2",dm);
	sf_putfloat(out,"o2",m0);

	/* for sffft3 */
	sf_putint(out,"m2",nh);
	sf_putfloat(out,"c2",dh);

	sf_putint(out,"n3",nk);
	sf_putfloat(out,"d3",dk);
	sf_putfloat(out,"o3",k0);

	if (NULL == sf_getstring("velocity")) {
	    /* file with velocity */
	    vel = NULL;
	} else {
	    vel = sf_input ("velocity");
	}
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
	    /* file with velocity */
	    vel = NULL;
	    if (!sf_getint("nz",&nz)) sf_error ("Need nz=");
	    /* Length of depth axis (for migration, if no velocity file) */
	    if (!sf_getfloat("dz",&dz)) sf_error ("Need dz=");
	    /* Sampling of depth axis (for migration, if no velocity file) */
	    if (!sf_getfloat("z0",&z0)) z0 = 0.;
	    /* Origin of depth axis (for migration, if no velocity file) */
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

    dw *= 2.0*SF_PI;
    w0 *= 2.0*SF_PI;
    dm *= 2.0*SF_PI;
    m0 *= 2.0*SF_PI;
    dk *= 2.0*SF_PI;
    k0 *= 2.0*SF_PI;

    vt = sf_floatalloc(nz);

    if (NULL == vel) {
	if (!sf_getfloat("vel",&v0)) sf_error ("Need vel=");
	/* Constant velocity (if no velocity file) */
	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	}
    } else {
	sf_floatread(vt,nz,vel);
    }

    /* allocate space */
    cp = sf_complexalloc2(nw,nm);
    cq = sf_complexalloc(nz);

    /* migrate each wavenumber */
    for (ik=0; ik<nk; ik++) {
	sf_warning("wavenumber %d of %d",ik+1,nk);
	
	k = ik*dk;

	if (inv) { /* modeling */
	    sf_complexread(cq,nz,in);
	} else {
	    for (iz=0; iz<nz; iz++) {
		cq[iz] = 0.0;
	    }
	    sf_complexread(cp[0],nw*nm,in);
	}

	for (im=0; im<nm; im++) {
	    jm = (im < nm/2)? im + nm/2: im - nm/2;
	    m = m0 + jm*dm;
      
	    dsr(inv,eps,k,m,nw,dw,w0,nz,dz,vt,cp[im],cq);
	}
 
	if (inv) {
	    sf_complexwrite(cp[0],nw*nm,out);
	} else {
	    sf_complexwrite(cq,nz,out);
	}
    }

    exit (0);
}

/* 	$Id$	 */
