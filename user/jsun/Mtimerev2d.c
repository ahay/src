/* 2-D correlative time reversal imaging of passive seismic data */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "absorb.c"

#define NOP 4 /* derivative operator half-size */
#define C0 -205.0f/72.0f
#define C1 +8.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 +8.0f/315.0f
#define C4 -1.0f/560.0f
#define Lap(a,ix,iz,sx,sz,v)  ( ( C4*(a[ix+4][iz  ] + a[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ] + a[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ] + a[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ] + a[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]) )*sx            +      \
                                ( C4*(a[ix  ][iz+4] + a[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3] + a[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2] + a[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1] + a[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]) )*sz )*v[ix][iz]

int main(int argc, char* argv[])
{
    bool verb,abc;                  /* execution flags */
    int ix, iz, it, ig;             /* index variables */
    int nt, wfnt, nx, nz, depth, nb, n2, snap, ngrp, counter;
    float ox, oz, dx, dz, dt, dt2, idz2, idx2, cb;
    int nxpad, nzpad;

    int *beg, *end;
    float **vvpad;
    float **dd, **vv, ***ww, ***wf;
    float **u0, **u1, **u2, **tmp;

    sf_file in, out, vel, wave;     /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);
    
    /* initialize OpenMP support */
#ifdef _OPENMP
    omp_init();
#endif

    if(!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
    if(!sf_getbool("abc",&abc)) abc = false;   /* absorbing boundary condition */
    if(!sf_getint("snap", &snap)) snap=0;      /* wavefield snapshot flag */
    if(!sf_getint("depth", &depth)) depth=0;   /* geophone depth */
    if(!sf_getint("ngrp", &ngrp)) ngrp=1;      /* number of groups */
	
    /* setup I/O files */
    in  = sf_input("in");
    out = sf_output("out");
    vel = sf_input("velocity");
    /* velocity model */
    
    /* Dimensions */
    if(!sf_histint  (vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histint  (vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "o1", &oz)) sf_error("No o1= in velocity");
    if(!sf_histfloat(vel, "o2", &ox)) sf_error("No o2= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");

    /* output files */
    if(!sf_histint(in, "n1", &nt)) sf_error("No n1= in data");
    if(!sf_histfloat(in, "d1", &dt)) sf_error("No d1= in data");
    if(!sf_histint(in, "n2", &n2) || n2!=nx) sf_error("Need n2=%d in data", nx);

    sf_putint   (out, "n1", nz);
    sf_putfloat (out, "o1", oz);
    sf_putfloat (out, "d1", dz);
    sf_putstring(out, "label1", "Depth");
    sf_putstring(out, "unit1" , "km");
    sf_putint   (out, "n2", nx);
    sf_putfloat (out, "o2", ox);
    sf_putfloat (out, "d2", dx);
    sf_putstring(out, "label2", "Distance");
    sf_putstring(out, "unit2" , "km");
    sf_putint   (out, "n3", nt);
    sf_putfloat (out, "d3", dt);
    sf_putfloat (out, "o3", 0.f);
    sf_putstring(out, "label3", "Time");
    sf_putstring(out, "unit3" , "s");

    /* dimension of padded boundary */
    if(!sf_getint("nb", &nb) || nb<NOP) nb = NOP;
    if(!sf_getfloat("cb", &cb)) cb = 0.0f;
    nxpad = nx+2*nb;
    nzpad = nz+2*nb;
    depth = depth+nb;

    /* set Laplacian coefficients */
    idz2 = 1.0f/(dz*dz);
    idx2 = 1.0f/(dx*dx);
    
    /* wavefield snapshot */
    if(snap){
        wfnt = 1+(nt-1)/snap;
        wave = sf_output("wave");
        
        sf_putint(wave, "n1", nz);
        sf_putint(wave, "d1", dz);
        sf_putint(wave, "o1", oz);
            
        sf_putint(wave, "n2", nx);
        sf_putint(wave, "d2", dx);
        sf_putint(wave, "o2", ox);
        
        sf_putint(wave, "n3", wfnt);
        sf_putfloat(wave, "d3", snap*dt);
        sf_putfloat(wave, "o3", 0.f);
    }
        
    /* allocate arrays */
    vv    = sf_floatalloc2(nz, nx);
    dd    = sf_floatalloc2(nt, nx);
    vvpad = sf_floatalloc2(nzpad, nxpad);
    u0    = sf_floatalloc2(nzpad, nxpad);
    u1    = sf_floatalloc2(nzpad, nxpad);
    u2    = sf_floatalloc2(nzpad, nxpad);
    ww = sf_floatalloc3(nz, nx, nt);
    beg=sf_intalloc(ngrp);
    end=sf_intalloc(ngrp);
    if (snap) wf = sf_floatalloc3(nz, nx, wfnt);
    else wf=NULL;
    
    /* read velocity */
    sf_floatread(vv[0], nz*nx, vel);
    
    /* pad boundary */
    dt2 = dt*dt;
    for     (ix=0; ix<nx; ix++)
        for (iz=0; iz<nz; iz++)
            vvpad[ix+nb][iz+nb] = vv[ix][iz]*vv[ix][iz]*dt2;
    for     (ix=0; ix<nxpad; ix++){
        for (iz=0; iz<nb;    iz++){
            vvpad[ix][      iz  ] = vvpad[ix][      nb  ];
            vvpad[ix][nzpad-iz-1] = vvpad[ix][nzpad-nb-1];
        }
    }
    for     (ix=0; ix<nb;    ix++){
        for (iz=0; iz<nzpad; iz++){
            vvpad[      ix  ][iz]=vvpad[      nb  ][iz];
            vvpad[nxpad-ix-1][iz]=vvpad[nxpad-nb-1][iz];
        }
    }

    /* absorbing boundary condition */
    if (abc) {
        if (verb) sf_warning("absorbing boundary condition");
        abc_init(nzpad,nxpad,nzpad,nxpad,nb,nb,nb,nb,cb,cb,cb,cb);
    }

    /* read data */
    sf_floatread(dd[0], nt*nx, in);

    /* set start and end index */
    counter = 0;
    for (ig=0; ig<ngrp; ig++) {
        beg[ig] = counter;
        counter += nx/ngrp;
        end[ig] = counter;
    }
    end[ngrp-1] = nx;
    if (verb) {
        for (ig=0; ig<ngrp; ig++) {
            sf_warning("beg[%d]=%d",ig,beg[ig]);
            sf_warning("end[%d]=%d",ig,end[ig]);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
    for         (it=0; it<nt; it++)
        for     (ix=0; ix<nx; ix++)
            for (iz=0; iz<nz; iz++)
                ww[it][ix][iz] = 1.0f;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
    for         (it=0; it<wfnt; it++)
        for     (ix=0; ix<nx; ix++)
            for (iz=0; iz<nz; iz++)
                wf[it][ix][iz] = 0.0f;

    for (ig=0; ig<ngrp; ig++) { /* loop over subgroups of receivers */

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<nxpad; ix++)
                for (iz=0; iz<nzpad; iz++)
                {
                    u0[ix][iz] = 0.0f;
                    u1[ix][iz] = 0.0f;
                    u2[ix][iz] = 0.0f;
                }

        for (it=nt-1; it>-1; it--){
            if (verb) sf_warning("Time reversal: %d/%d;", it, 0);

            /* time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<nxpad-NOP; ix++){
                for (iz=NOP; iz<nzpad-NOP; iz++){
                    u2[ix][iz] = Lap (u1,ix,iz,idx2,idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;
            if (abc) abc_apply(u1[0]);
            if (abc) abc_apply(u0[0]);

            /* inject data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix=nb+beg[ig]; ix<nb+end[ig]; ix++)
                u1[ix][depth] += dd[ix-nb][it];

            /* image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<nx; ix++)
                for (iz=0; iz<nz; iz++)
                    ww[it][ix][iz] *= u1[ix+nb][iz+nb];

            if (snap && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
                for     (ix=0; ix<nx; ix++)
                    for (iz=0; iz<nz; iz++)
                        wf[it/snap][ix][iz] += u1[ix+nb][iz+nb];
            }

        } /* it loop */
        if (verb) sf_warning(".");

    } /* ig loop */

    /* output source */
    sf_floatwrite(ww[0][0], nz*nx*nt, out);
    if (snap) sf_floatwrite(wf[0][0], nz*nx*wfnt, wave);

    free(**ww); free(*ww); free(ww);
    if (snap) { free(**wf); free(*wf); free(wf); }

    if (abc) abc_close();
    free(*vvpad); free(vvpad); free(*vv); free(vv);
    free(*dd); free(dd); 
    free(*u0); free(u0); free(*u1); free(u1); free(*u2); free(u2);
    exit (0);
}

