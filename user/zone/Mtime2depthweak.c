/* Time-to-depth conversion in media with weak lateral variations 2D (Sripanich and Fomel, 2017). */

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

int main (int argc, char *argv[])
{
    int nz, nx, iz, ix, izz, subsample, k, nderiv, nsmooth, smoothlen;
    float  z0, dz, x0, dx, c, deldistx, deltimex, refvelmid, velmid, dvdx0mid, dvdt0mid, deldistmid,deltimemid;
    float **vel, **dvdx0, **dvdt0, **refvel, **deltime, **deldist, **delv, *currentx, *currentt, *currentv, *dcurrentx, *dcurrentt, *dcurrentv;
    sf_triangle smooth;
    sf_file in, dveldx0, dveldt0, refvelocity, outdeltime, outdeldist, outdelv;

    sf_init(argc, argv);
    in          = sf_input("in"); /* velocity squared from preliminary dix inversion */
    dveldx0     = sf_input("dvdx0"); /* from derivative of velocity squared in time-domain x0 + preliminary dix inversion */
    dveldt0     = sf_input("dvdt0"); /* from derivative in of velocity squared time-domian t0 + from preliminary dix inversion */
    refvelocity = sf_input("refvelocity"); /* ref velocity squared from preliminary dix inversion */
    
    outdeltime = sf_output("outdt0"); /* dt0*/
    outdeldist = sf_output("outdx0"); /* dx0*/
    outdelv    = sf_output("outdv");  /* dv*/

    if (!sf_histint (in,"n1",&nz)) sf_error ("No n1= in input");
    if (!sf_histfloat (in,"d1",&dz)) sf_error ("No d1= in input");
    if (!sf_histfloat (in,"o1",&z0)) z0 = 0.;
    
    if (!sf_histint (in,"n2",&nx)) sf_error ("No n2= in input");
    if (!sf_histfloat (in,"d2",&dx)) sf_error ("No d2= in input");
    if (!sf_histfloat (in,"o2",&x0)) x0 = 0.0 ;

    if (!sf_getint ("zsubsample",&subsample)) subsample = 100;
    /* Additional subsampling in depth for stability */

    vel = sf_floatalloc2 (nz,nx);
    dvdx0 = sf_floatalloc2 (nz,nx);
    dvdt0 = sf_floatalloc2 (nz,nx);
    
    refvel = sf_floatalloc2 (nz,nx);
    deltime = sf_floatalloc2(nz,nx);
    deldist = sf_floatalloc2(nz,nx);
    delv = sf_floatalloc2(nz,nx);
    
    currentx = sf_floatalloc(nx);
    currentt = sf_floatalloc(nx);
    currentv = sf_floatalloc(nz);
    dcurrentx = sf_floatalloc(nx);
    dcurrentt = sf_floatalloc(nx);
    dcurrentv = sf_floatalloc(nz);
    

	/* Read input*/
	sf_floatread (vel[0],nz*nx,in);
	sf_floatread (dvdx0[0],nz*nx,dveldx0);
	sf_floatread (dvdt0[0],nz*nx,dveldt0);
	sf_floatread (refvel[0],nz*nx,refvelocity);
	
	/* Setup derivative filter*/
	
    if (!sf_getint("nderiv",&nderiv)) nderiv=10;
    /* Derivative filter order */
    if (!sf_getfloat("refderiv",&c)) c=1.;
    /* Deriveative filter reference (0.5 < ref <= 1) */
	
	sf_deriv_init(nx, nderiv, c); 
	
	if (!sf_getint("smoothlen",&smoothlen)) smoothlen=nx/20;
    /* Smoothing filter length */
    if (!sf_getint("nsmooth",&nsmooth)) nsmooth=10;
    /* Smoothing repeat */
	smooth = sf_triangle_init(smoothlen,nx,false);
	
	
	/* Weak variation contribution (stepping in depth not optimal efficiency in the current array allocation)*/
	for (iz = 0; iz < nz-1; iz++) {
		for (k = 0; k < nx; k++) {
				currentx[k] = deldist[k][iz];
				currentt[k] = deltime[k][iz];
		}
		
		/* Find the x derivative and smooth*/
		sf_deriv(currentx,dcurrentx);
		sf_deriv(currentt,dcurrentt);
		
		for (k = 0 ; k < nsmooth ; k++ ) sf_smooth2(smooth,0,1,false,dcurrentx);
		for (k = 0 ; k < nsmooth ; k++ ) sf_smooth2(smooth,0,1,false,dcurrentt);
		
		for (ix = 0; ix < nx; ix++) {
			deltimex = dcurrentt[ix]/dx;
			deldistx = dcurrentx[ix]/dx;
			deltimemid = deltime[ix][iz];
			deldistmid = deldist[ix][iz];
		// Finer depth stepping for stability
			for (izz = 0; izz < subsample; izz++) {
				refvelmid = refvel[ix][iz] + izz*(refvel[ix][iz+1] - refvel[ix][iz])/subsample;
				velmid    = vel[ix][iz] + izz*(vel[ix][iz+1] - vel[ix][iz])/subsample;
				dvdx0mid  = dvdx0[ix][iz] + izz*(dvdx0[ix][iz+1] - dvdx0[ix][iz])/subsample;
				dvdt0mid  = dvdt0[ix][iz] + izz*(dvdt0[ix][iz+1] - dvdt0[ix][iz])/subsample;
				
				deltimemid = deltimemid + (dz/subsample)*( deldistx/sqrt(refvelmid) - (1/(2*pow(refvelmid,1.5)))*(dvdt0mid*deltimemid + dvdx0mid*deldistmid + velmid - refvelmid) );
				deldistmid = deldistmid + (dz/subsample)*( -sqrt(refvelmid)*deltimex );
			}
			deltime[ix][iz+1] = deltimemid;
			deldist[ix][iz+1] = deldistmid;
		}
		sf_warning("Depth step: %d of %d;",iz+1,nz);
	}
	sf_warning(".");
	
	// Compute dv2
	sf_deriv_init(nz, nderiv, c); 
	for (ix = 0; ix < nx; ix++) {
		for (iz = 0; iz < nz; iz++) {
			currentv[iz] = deltime[ix][iz];
		}
		sf_deriv(currentv,dcurrentv);
		for (iz = 0; iz < nz; iz++) {
			delv[ix][iz] = -2*pow(refvel[ix][iz],1.5)*dcurrentv[iz]/dz;
		}
	}

	sf_floatwrite (deltime[0],nz*nx,outdeltime);
	sf_floatwrite (deldist[0],nz*nx,outdeldist);
	sf_floatwrite (delv[0],nz*nx,outdelv);

    free(vel); free(dvdx0); free(dvdt0);
    free(refvel); free(deltime); free(deldist);
    free(delv);
    
    exit (0);
}

/* 	$Id$	 */
