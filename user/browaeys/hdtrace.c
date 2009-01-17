/* Multiple arrivals by marching down/up using Hamiltonian dynamics. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include <math.h>
#include <rsf.h>

#include "hdtrace.h"
#include "symplectrace.h"

#ifdef _it_does_not_compile

include "grad2fill.h"

#ifndef _hdtrace_h

#define NS 4
/*^*/

#endif

static int nx, nz, np, npx;
static float dx, dz, dp, x0, z0, p0, xm, zm, pm;
static sf_eno2 cvel, fslc[NS];
static float **slice;
static bool *known;
static pqv hvec;

static const float eps = 1.e-5;
static void psnap (float* p, float* q, int* iq);


void hdtrace_init (int order        /* interpolation order for velocity */, 
		   int iorder       /* interpolation order for values */,
		   int nx1          /* lateral samples */, 
		   int nz1          /* depth samples */, 
		   int np1          /* horizontal slowness samples */, 
		   float dx1        /* lateral sampling */, 
		   float dz1        /* depth sampling */, 
		   float dp1        /* horizontal slowness sampling */,
		   float x01        /* lateral origin */, 
		   float z01        /* depth origin */, 
		   float p01        /* horizontal slowness origin */,
                   float** vel      /* slowness [nx][nz] */,
                   float** slice_in /* depth slice [nx][np][NS] */)
/*< Initialize >*/
{
    int is;

    slice = slice_in;

    nx = nx1; 
    nz = nz1; 
    np = np1;

    dx = dx1;
    dz = dz1;
    dp = dp1;

    x0 = x01; 
    z0 = z01;
    p0 = p01;

    npx = np*nx; 

    xm = x0 + (nx-1)*dx;
    zm = z0 + (nz-1)*dz;
    pm = p0 + (np-1)*dp;

    cvel = sf_eno2_init (order,nz,nx);
    sf_eno2_set (cvel,vel);

    known = sf_boolalloc (npx);

    grad2fill_init (np,nx);

    for (is=0; is < NS; is++) {
        fslc[is] = sf_eno2_init (iorder,np,nx);
    }

    nc4_init();
    hvec = (pqv)malloc(sizeof(pqv));
}

void hdtrace_step (int kz, int up, float **slow) 
/*< Step in depth >*/
{
    int incell, is, nk, kx, kp, k, ix, iz, step;
    float z, x, p[2], t, ds;
    float ssi, gi[2], ss, g[2];
    float fx, fz, fp;
    bool onx, onz;

    /* grid dimension restrictions */
    /* dx.dp > 1/(4.pi.f) */
    /* dz.dp > 1/(4.pi.f) */

    if (up == -1) sf_warning("upgoing rays");
    if (up ==  1) sf_warning("downgoing rays");

    for (is=0; is < NS; is++) {    
        sf_eno2_set1 (fslc[is],slice[is]);
    }

    nk = 0; /* number of known */

    z = z0 + kz*dz;

    for (kx = 0; kx < nx; kx++) {

	x = x0 + kx*dx;

        sf_eno2_apply(cvel,kz,kx,0.,0.,&ssi,gi,BOTH);
        gi[1] /= dx;
        gi[0] /= dz;

        for (kp = 0; kp < np; kp++) {

            k = kp + kx*np;

           /* initial dimensionless horizontal slowness */
	    p[1] = p0 + kp*dp;

            /* initial dimensionless vertical one-way slowness */
	    p[0] = up*sqrt(1.-p[1]*p[1]);

            ss = ssi;
            g[0] = gi[0];
            g[1] = gi[1];

            ix = kx; 
	    onx = true;

            iz = kz; 
	    onz = true;

	    incell = 1;

            /* decide if we are out already */
            if ((iz == 0    && p[0] < 0.) ||
                (iz == nz-1 && p[0] > 0.) ||
                (ix == 0    && p[1] < 0.) ||
                (ix == nx-1 && p[1] > 0.)) {
                slice[0][k] = 0.;
                slice[1][k] = x0 + ix*dx;
                slice[2][k] = z0 + iz*dz;
                slice[3][k] = p[1];
                known[k] = true;
                nk++;
                continue;
            } else {
                known[k] = false;
            }

            /* Hamiltonian vector */
	    hvec_init(hvec,0.,x,z,ss0*p[1],ss0*p[0]);

            /* predict sigma step size to the next cell */
	    ds = nc4_cellstep(hvec,slow,nx,nz,dx,dz,x0,z0,dp,dp);

            /* Loop until previous depth level or exiting from phase space grid */
	    while (incell == 1) {

                /* symplectic cell step and traveltime integration */
		nc4_sigmastep(hvec,fabs(ds),&ss,slow,nx,nz,dx,dz,x0,z0);

		value_exitlevel(hvec,ss,&step,&fx,&fz,&fp,&t);

                /* exit at previous/next depth level */
		if (step == 1) {
		    onz = true;
                    onx = sf_cell_snap (&fx,&ix,eps);
		    incell = 0;
		    break;
		}

                /* exit from spatial grid sides */
		if (step == 2) {
		    onx = true;
		    onz = sf_cell_snap (&z,&iz,eps);
                    if ((xm-fx) <= dx && fp > 0.) {
                        /* exit from the right side */

			enogrid_apply(slice,nx-1,0.,p2,f);
			f[2] = z1;
			incell = 0;
			break;
		    }
		    else if (fx <= dx && fp < 0.) {
                        /* exit from the left side */
			enogrid_apply(slice,0,0.,p2,f);
			f[2] = z1;
			incell = 0;
			break;
		    }	        
		}

      ix < 0 || ix > nx-1 ||
                    (onx && ix == 0 && p[1] < 0.) ||
                    (ix == nx-1 && (!onx || p[1] > 0.))) {
                    slice[0][k] = t;
                    slice[1][k] = x0+(x+ix)*dx;
                    slice[2][k] = z0+(z+iz)*dz;
                    slice[3][k] = sf_cell_p2a(p)*180./SF_PI;
                    known[k] = true;
                    nk++;
                    break;



                /* update for next sigma step */
		ds = nc4_cellstep(hvec,slow,nx,nz,dx,dz,x0,z0,dp,dp);

                /* exit from slowness grid sides (overturning ray) */
		if (step == 3) {
		    if ((pm-p2) <= dp && ds > 0.) {
			fx = (x2-x0)/dx;
			ix = snap(&x2,nx);
			enogrid_apply(slice,ix,fx,p2,f);
			f[2] = z1;
			incell = 0;
			break;
		    }
                    else if (p2 <= dp && ds < 0.) {
			fx = (x2-x0)/dx;
			ix = snap(&x2,nx);
			enogrid_apply(slice,ix,fx,p2,f);
			f[2] = z1;
			incell = 0;
			break;
		    }
		}
	    }


	    slice[0][k] = t;
	    slice[1][k] = x0+(x+ix)*dx;
	    slice[2][k] = z0+(z+iz)*dz;
	    slice[3][k] = p[1];
	    known[k] = true;
	    nk++;

	} /* kp */
    } /* kx */

    if (nk < npx) {
        fprintf(stderr,"known=%d (%d)\n",nk,npx);
        nk = SF_MIN(SF_MIN(np,nx),npx-nk);
        for (is=0; is < NS; is++) {
            grad2fill (nk, slice[is], known):
        }
    }

}





static void psnap (float* p, float* q, int* iq)
/* snap to the grid (non dimensional horizontal slowness) */
{
    int ip;
    float fp2, fp;

    fp = (p-p0)/dp;
    ip = floor (fp); 
    fp2 = fp-ip;
    sf_cell_snap (&fp2,&ip,eps);

    if (ip < 0) {
        ip = 0; 
	fp2 = 0.;
	fp = p0;
    } else if (ip >= np || (ip == np-1 && fp2 > 0.)) {
        ip = np-1;
	fp2 = 0.;
	fp = f0+(np-1)*dp;
    }

    p[1] = fp;              /* px is  sin(a) */
    p[0] = -sqrt(1.-fp*fp); /* pz is -cos(a) */

    *q = fp2;
    *iq = ip;
}

#endif
