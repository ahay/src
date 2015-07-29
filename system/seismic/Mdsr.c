/* Prestack 2-D VTI v(z) modeling/migration by DSR with angle gathers. */
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
    int nt, nt2;               /* number of time samples */
    int nz;		/* number of migrated time samples */
    int nw;               /* number of frequencies */
    int nh;               /* number of offsets */
    int ik,im,iz,ia;      /* loop counters 	*/
    int nk,nm;		/* number of wave numbers */	
    int na;             /* number of angles */

    float dt;             /* time sampling interval */
    float dz;		/* migrated time sampling interval */
    float dh;             /* offset sampling interval */
    float dk,dm;	        /* wave number sampling interval */
    float da;              /* angle sampling interval */
    float k,m;            /* wave number                  */
    float k0=0.;       /* first offset wavenumber */
    float *vt, v0;	/* velocity v(t)		*/
    float *vz, vz0;	/* vertical velocity v(t)		*/
    float *n, n0;	/* eta		*/

    float *p, **q;	/* data, image */

    bool inv;              /* modeling or migration        */
    bool depth;           /* time or depth migration      */
    float eps;            /* dip filter constant          */   
    
    const char *rule;   /* phase-shuft interpolation rule */
    const char *arule;    /* angle gather rule */
            
    sf_file vel, velz=NULL, eta=NULL, in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; If n, migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* Stabilization parameter */

    if (!sf_getbool("depth",&depth)) depth = false;
    /* if true, depth migration */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("na",&na)) na=1;
    /* number of angles */

    if (!sf_getfloat("da",&da)) da=90.;
    /* angle sampling (in degrees) */

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");

	if (!sf_histint(in,"n2",&nk)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dk)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&k0)) k0=0.;

	if (!sf_getint("nt",&nt)) {
	    /* Length of time axis (for modeling) */
	    if (depth) {
		sf_error ("Need nt=");
	    } else {
		nt = nz;
	    }
	} else {
	    sf_putint(out,"n1",nt);
	}

	if (!sf_getfloat("dt",&dt)) {
	    /* Sampling of time axis (for modeling) */
	    if (depth) {
		sf_error ("Need dt=");
	    } else {
		dt = dz;
	    }
	} else {
	    sf_putfloat(out,"d1",dt);
	}

	sf_putfloat(out,"o1",0.);

	/*
	  if (!sf_getfloat("t0",&t0)) t0 = 0.;
	  sf_putint(out,"nt",nt);
	  sf_putfloat(out,"t0",t0);

	  nw = nt;
	  dw = 1./(nw*dt);
	  w0 = -0.5/dt;
	  
	  sf_putint(out,"n1",nw);
	  sf_putfloat(out,"d1",dw);
	  sf_putfloat(out,"o1",w0);
	*/	  

	if (!sf_getint("nh",&nh)) sf_error ("Need nh=");
	/* Number of offsets (for modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error ("Need dh=");
	/* Offset sampling (for modeling) */

	/* determine wavenumber sampling, pad by 2 */
	nm = nh*2;
	nm *= 2;
	dm = 1./(2*(nm-1)*dh);

	sf_putint(out,"n2",nm);
	sf_putfloat(out,"d2",dm);
	sf_putfloat(out,"o2",0.);

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
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1= in input");

	if (!sf_histint(in,"n2",&nm)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dm)) sf_error ("No d2= in input");
	
	if (!sf_histint(in,"n3",&nk)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dk)) sf_error ("No d3= in input");
	if (!sf_histfloat(in,"o3",&k0)) sf_error ("No d3= in input");

	if (NULL == sf_getstring("velocity")) {
	    /* file with velocity */
	    vel = NULL;
	    velz = NULL;
	    eta = NULL;
	    if (!sf_getint("nz",&nz)) nz=nt;
	    /* Length of depth axis (for migration, if no velocity file) */
	    if (!sf_getfloat("dz",&dz)) dz=dt;
	    /* Sampling of depth axis (for migration, if no velocity file) */
	} else {
	    vel = sf_input ("velocity");
	    if (!sf_histint(vel,"n1",&nz)) sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) sf_error ("No d1= in velocity");

	    if (NULL == sf_getstring("velz")) {
		velz = NULL;
		eta = NULL;
	    } else {
		velz = sf_input("velz");
		if (NULL == sf_getstring("eta")) sf_error("Need eta=");
		eta = sf_input("eta");
	    }

	    
	}

	sf_putint(out,"n1",na);
	sf_putfloat(out,"d1",da);
	sf_putfloat(out,"o1",0.);

	sf_putint(out,"n2",nz);
	sf_putfloat(out,"d2",dz);
	sf_putfloat(out,"o2",0.);
    }

    da *= SF_PI/180.; /* convert to radians */

    dm *= 2.0*SF_PI;
    dk *= 2.0*SF_PI;
    k0 *= 2.0*SF_PI;

    vt = sf_floatalloc(nz);
    vz = sf_floatalloc(nz);
    n = sf_floatalloc(nz);

    if (NULL == (rule = sf_getstring("rule"))) rule="simple";
    /* phase-shift interpolation rule (simple, midpoint, linear, anisotropic, dti) */

    if (NULL == (arule = sf_getstring("arule"))) arule="isotropic";
    /* angle gather rule */

    if (NULL == vel) {
	if (!sf_getfloat("vel",&v0)) sf_error ("Need vel=");
	/* Constant velocity (if no velocity file) */

	if (!sf_getfloat("vz",&vz0)) vz0=v0;
	/* Constant vertical velocity (if no velocity file) */

	if (!sf_getfloat("n",&n0)) n0=0.0;
	/* Constant eta (if no velocity file) */

	for (iz=0; iz < nz; iz++) {
	    vt[iz] = v0;
	    vz[iz] = vz0;
	    n[iz] = n0;
	}
    } else {
	sf_floatread(vt,nz,vel);
	sf_fileclose(vel);

	if ('a' == rule[0] || 'a' == arule[0] || 'd' == rule[0] || 'd' == arule[0] ) {
	    if (NULL == velz || NULL == eta) sf_error("Need velz= and eta=");
	    sf_floatread(vz,nz,velz);
	    sf_floatread(n,nz,eta);

	    sf_fileclose(velz);
	    sf_fileclose(eta);
	}  else {
	    for (iz=0; iz < nz; iz++) {
		vz[iz] = vt[iz];
		n[iz] = 0.0;
	    }
	}
    }

    /* vt -> vt^2 */ 
    for (iz=0; iz < nz; iz++) {
	vt[iz] *= vt[iz];
	vz[iz] *= vz[iz];
	if (depth) { /* convert to slowness */
	    vt[iz] = 1./vt[iz];
	    vz[iz] = 1./vz[iz];
	}
    }

    /* determine frequency sampling */    
    nt2 = 2*kiss_fft_next_fast_size((nt+1)/2);

    if (!sf_getint("nw",&nw)) nw = nt2/2+1;
    /* Maximum number of frequencies */

    /* allocate space */
    p = sf_floatalloc(nt2);
    q = sf_floatalloc2(na,nz);
    
    dsr_init(eps, nw, nt2, dt, nz, dz, vt, vz, n, depth, rule[0], na, da);
    
    /* migrate each wavenumber */
    for (ik=0; ik<nk; ik++) { /* midpoint wavenumber */
	sf_warning("wavenumber %d of %d",ik+1,nk);
	
	k = k0+ik*dk;

	if (inv) {
	    sf_floatread(q[0],nz*na,in);
	} else {
	    /* initialize image to zero */
	    for (iz=0; iz < nz; iz++) {
		for (ia=0; ia < na; ia++) {
		    q[iz][ia] = 0.;
		}
	    }
	}

	for (im=0; im<nm; im++) { /* half-offset wavenumber */
	    m = im*dm;      

	    if (!inv) {
		sf_floatread(p,nt,in);
		if (nt != nt2) p[nt] = 0.;
	    }

	    dsr(arule[0],inv,k,m,p,q);

	    if (inv) sf_floatwrite(p,nt,out);
	}
 
	/* output image */
	if (!inv) sf_floatwrite(q[0],nz*na,out);
    }


    exit (0);
}

/* 	$Id: Mdsr.c 7590 2011-08-13 01:33:29Z ivlad $	 */
