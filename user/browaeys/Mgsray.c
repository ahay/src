/* Gauss Seidel iterative solver for phase space escape positions, angle and traveltime */
/*
   Downwind first-order scheme with angular grid ]-180,180]
   Escape angle from vertical for slowness pointing upward

   Equation for escape positions or angle is (with appropriate BC)
   -S.sin(a).dX/dx - S.cos(a).dX/dz - [cos(a).Sx - sin(a).Sz].dX/da = 0

   Equation for escape traveltime T is (with appropriate BC)
   -S.sin(a).dT/dx - S.cos(a).dT/dz - [cos(a).Sx - sin(a).Sz].dT/da = -S*S
*/
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include <assert.h>

#include "gsray.h"
#include "eno3.h"

int main(int argc, char* argv[])
{
    int nz,nx,na;
    float oz,ox,oa;
    float dz,dx,da;

    int iq;                        /* escape variables switch */  

    int ix;                        /* grid points in x */
    int iz;                        /* grid points in z */
    int ia;                        /* grid points in a */

    int iter, niter, liter;        /* number of iterations */
    int order;                     /* order of upwind */

    int ix0, ix1, ixs;
    int iz0, iz1, izs;
    int ia0, ia1, ias;

    float a;                 /* angle */

    float ***t;              /* escape variable */
    float **s,**sx,**sz;     /* slowness, gradients */
    float **ds,**dsx,**dsz;  /* delta slowness, gradients */
    float ***rhs;            /* extra RHS for linearization */
    eno3 dt0;                /* interpolation for T0 in angle in extra RHS */
    int step, nsteps = 1;
    enum { INIT=0, LINEAR=1 }; /* Possible calculation steps, 0 - nonlinear(initial)
                                                              1 - further linearization with dslow */

    float ss,ssx,ssz;
    float dss,dssx,dssz;

    float trhs,lrhs,arhs;
    float grad[3];

    float cs, sn, zn, zf, xn, xf, an, af;
    float new_val;
    float cvt,tol,cvt_old = SF_HUGE;
    bool verb, sph;

    sf_file in,out,slow,slowz,slowx,dslow,dtout=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("iq",&iq)) sf_error("Need iq=");
    /* switch for escape variable 0=x, 1=a, 2=t, 3=z, 4=l, 5=i */

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&oz)) sf_error("No o1=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");

    if (!sf_histint(in,"n3",&na)) sf_error("No n3=");
    if (!sf_histfloat(in,"d3",&da)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&oa)) sf_error("No o3=");
    /* angle in degrees */

    if (!sf_getint("niter",&niter)) niter=50;
    /* number of Gauss-Seidel iterations */
    if (!sf_getint("liter",&liter)) liter=0;
    /* number of first iterations with low-order scheme */

    if (!sf_getfloat("tol",&tol)) tol=0.000002*nx*nz;
    /* accuracy tolerance */

    if (!sf_getint("order",&order)) order=1;
    /* order of upwind */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool ("sph", &sph)) sph = false;
    /* true - half-sphere, false - flat B.C. on left/right */

    /* memory allocations */

    /* Solution vector */
    t = sf_floatalloc3(nz,nx,na);

    /* read input escape variable - initial guess */
    sf_floatread(t[0][0],nz*nx*na,in);

    /* read auxiliary slowness file */
    slow = sf_input("slow");
    s = sf_floatalloc2(nz,nx);
    sf_floatread(s[0],nz*nx,slow);

    /* read auxiliary slowness z-gradient file */
    slowz = sf_input("slowz");
    sz = sf_floatalloc2(nz,nx);
    sf_floatread(sz[0],nz*nx,slowz);

    /* read auxiliary slowness x-gradient file */
    slowx = sf_input("slowx");
    sx = sf_floatalloc2(nz,nx);
    sf_floatread(sx[0],nz*nx,slowx);

    rhs=NULL;
    ds=dsx=dsz=NULL;
    if (sf_getstring("dslow") && sf_getstring("dtout")) {
        dslow = sf_input("dslow");
        ds = sf_floatalloc2(nz,nx);
        dsx = sf_floatalloc2(nz,nx);
        dsz = sf_floatalloc2(nz,nx);
        sf_floatread(ds[0],nz*nx,dslow);
        sf_floatread(dsz[0],nz*nx,dslow);
        sf_floatread(dsx[0],nz*nx,dslow);
        nsteps = 2;
        rhs = sf_floatalloc3(nz,nx,na);
        dtout = sf_output("dtout");
    }

    /* convert to radians */
    oa *= SF_PI/180.;
    da *= SF_PI/180.;

    gsray_init(nz,nx,na,oz,ox,oa,dz,dx,da);
    gsray_order(order,liter);

    for (step = INIT; step < nsteps; step++) {
    ia0=0; ia1=na; ias=1;
    for (iter=0; iter < niter; iter++) {
        cvt = 0.0; /* X_next - X_prev */

        /* alternate updating direction on angle grid*/
        if (iter % 2) {
            ia0 = 0; ia1 = na; ias = 1;
        } else {
            ia0 = na-1; ia1 = -1; ias = -1;
        }
        /* Gauss-Seidel iteration on angle */
        for (ia = ia0; ia != ia1; ia += ias) {
            a = oa + ia*da;
            cs = cosf(a);
            sn = sinf(a);
            assert (fabs(sn) > 1e-6 && fabs(cs) > 1e-6);

            if (1e-6 < -sn) {
                ix0=nx-2; ix1=-1; ixs=-1;
                for (iz = 0; iz < nz; iz++) {
                    gs_init_bc(t,s,sph,iz,nx-1,ia,INIT == step ? iq : 6);
                }
            } else {
                ix0=1; ix1=nx; ixs=1;
                for (iz = 0; iz < nz; iz++) {
                    gs_init_bc(t,s,sph,iz,0,ia,INIT == step ? iq : 6);
                }
            }
            if (1e-6 < -cs) {
                iz0=nz-2; iz1=-1; izs=-1;
                for (ix = 0; ix < nx; ix++) {
                    gs_init_bc(t,s,sph,nz-1,ix,ia,INIT == step ? iq : 6);
                }
            } else {
                iz0=1; iz1=nz; izs=1;
                for (ix = 0; ix < nx; ix++) {
                    gs_init_bc(t,s,sph,0,ix,ia,INIT == step ? iq : 6);
                }
            }
            /* loop over grid Z.X */
            for (ix = ix0; ix != ix1; ix += ixs) {
                for (iz = iz0; iz != iz1; iz += izs) {
                    ss = s[ix][iz];
                    ssx = sx[ix][iz];
                    ssz = sz[ix][iz];

                    if (LINEAR == step) {
                        dss = ds[ix][iz];
                        trhs = rhs[ia][ix][iz];
                        lrhs = dss - trhs;
                        arhs = -trhs;
                        trhs = 2.0*ss*dss - trhs;
                    } else {
                        trhs = ss*ss; /* RHS for escape time */
                        lrhs = ss; /* RHS for ray length */
                        arhs = 0.; /* RHS for escape x,z,a */
                    }
                    /* Gauss-Seidel update */
                    if (order > 3)
                        gs_get_face_cf (s, sz, sx, iz, ix, ia, &zn, &zf, &xn, &xf, &an, &af);
                    new_val = gs_update (t,cs*ss, sn*ss, (cs*ssx - sn*ssz), iz,ix,ia, trhs, lrhs, arhs, iq,
                                         iter, zn, zf, xn, xf, an, af);

                    cvt += fabsf(t[ia][ix][iz]-new_val);
                    t[ia][ix][iz] = new_val;
                } /* ix */
            } /* iz */
        } /* ia */
        if (verb) sf_warning("Iter=%d, Norm L1=%g",iter,cvt/(nx*nz*na));
        /* tol is tolerance for convergence */
        if ((cvt_old - cvt) >= 0 && (cvt_old - cvt) < tol && iter >= liter) break;
        cvt_old = cvt;
    } /* end G-S iterations */
    if (INIT == step) {
        sf_floatwrite(t[0][0],nz*nx*na,out);
    }
    /* Prepare RHS for linearization */
    if (nsteps > 1 && INIT == step) {
        if (verb) sf_warning("Preparing RHS for linearization");
        dt0 = eno3_init(3,nz,nx,na);
        eno3_set(dt0,t);
        for (ia = 0; ia < na; ia++) {
           a = oa + ia*da;
           cs = cosf(a);
           sn = sinf(a);
           for (ix = 0; ix < nx; ix++) {
                for (iz = 0; iz < nz; iz++) {
                    eno3_apply(dt0,iz,ix,ia,0.,0.,0.,&trhs,grad,BOTH);
                    dss = ds[ix][iz];
                    dssx = dsx[ix][iz];
                    dssz = dsz[ix][iz];
                    rhs[ia][ix][iz] = (grad[0]/dz*cs + grad[1]/dx*sn)*dss + grad[2]/da*(dssx*cs - dssz*sn);
                }
            }
        }
        eno3_close(dt0);
        memset(&t[0][0][0],0,na*nx*nz*sizeof(float));
        if (verb) sf_warning("Running linearization");
    }
    /* output */
    if (LINEAR == step)
        sf_floatwrite(t[0][0],nz*nx*na,dtout);
    } /* end steps */
    
    exit(0);
}


