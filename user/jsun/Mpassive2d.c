/* 2-D passive seismic RTM and its adjoint */
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
#define LapT(a,ix,iz,sx,sz,v) ( ( C4*(a[ix+4][iz  ]*v[ix+4][iz  ] + a[ix-4][iz  ]*v[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ]*v[ix+3][iz  ] + a[ix-3][iz  ]*v[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ]*v[ix+2][iz  ] + a[ix-2][iz  ]*v[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ]*v[ix+1][iz  ] + a[ix-1][iz  ]*v[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sx                          +      \
                                ( C4*(a[ix  ][iz+4]*v[ix  ][iz+4] + a[ix  ][iz-4]*v[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3]*v[ix  ][iz+3] + a[ix  ][iz-3]*v[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2]*v[ix  ][iz+2] + a[ix  ][iz-2]*v[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1]*v[ix  ][iz+1] + a[ix  ][iz-1]*v[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sz )

int main(int argc, char* argv[])
{
    bool verb,pas,adj,abc;          /* execution flags */
    int ix, iz, it;                 /* index variables */
    int nt, nx, nz, depth, nzxpad, nb, n2, snap;
    float ox, oz, dx, dz, dt, dt2, idz2, idx2, cb;

    int nxpad, nzpad;
    float **vvpad;
    
    float **dd, **mm, **vv, ***ww;
    float **u0, **u1, **u2, **tmp;  /* temporary arrays */

    sf_file in, out, vel, wave=NULL;     /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);
    
    /* initialize OpenMP support */
#ifdef _OPENMP
    omp_init();
#endif

    if(!sf_getbool("verb", &verb)) verb=false;
    /* verbosity flag */
    if(!sf_getbool("adj", &adj)) adj=false;
    /* adjoint flag, 0: modeling, 1: migration */
    if(!sf_getbool("pas", &pas)) pas=false;
    /* passive flag, 0: exploding reflector rtm, 1: passive seismic imaging */
    if(!sf_getbool("abc",&abc)) abc = false;
    /* absorbing boundary condition */
    if(!sf_getint("snap", &snap)) snap=0;
    /* wavefield snapshot flag */
    if(!sf_getint("depth", &depth)) depth=0;
    /* surface */
	
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

    if(adj){ /* migration */
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
        if (pas) {
            sf_putint   (out, "n3", nt);
            sf_putfloat (out, "d3", dt);
            sf_putfloat (out, "o3", 0.0f);
            sf_putstring(out, "label3", "Time");
            sf_putstring(out, "unit3" , "s");
        }
    }else{ /* modeling */
        if(!sf_getint("nt", &nt)) sf_error("Need nt=");
        if(!sf_getfloat("dt", &dt)) sf_error("Need dt=");
        
        sf_putint   (out, "n1", nt);
        sf_putfloat (out, "d1", dt);
        sf_putfloat (out, "o1", 0.0);
        sf_putstring(out, "label1", "Time");
        sf_putstring(out, "unit1" , "s");
        sf_putint   (out, "n2", nx);
        sf_putfloat (out, "o2", ox);
        sf_putfloat (out, "d2", dx);
        sf_putstring(out, "label2", "Distance");
        sf_putstring(out, "unit2" , "km");
        if (pas) {
            sf_putint   (out, "n3", 1);
        }
    }
    
    /* dimension of padded boundary */
    if(!sf_getint("nb", &nb) || nb<NOP) nb = NOP;
    if(!sf_getfloat("cb", &cb)) cb = 0.0f;
    nxpad = nx+2*nb;
    nzpad = nz+2*nb;
    nzxpad = nzpad*nxpad;
    depth = depth+nb;

    /* set Laplacian coefficients */
    idz2 = 1.0f/(dz*dz);
    idx2 = 1.0f/(dx*dx);
    
    /* wavefield snapshot */
    if(snap){
        wave = sf_output("wave");
        
        sf_putint(wave, "n1", nzpad);
        sf_putfloat(wave, "d1", dz);
        sf_putfloat(wave, "o1", oz-nb*dz);
            
        sf_putint(wave, "n2", nxpad);
        sf_putfloat(wave, "d2", dx);
        sf_putfloat(wave, "o2", ox-nb*dx);
        
        sf_putint(wave, "n3", 1+(nt-1)/snap);
        if(adj){
            sf_putfloat(wave, "d3", -snap*dt);
            sf_putfloat(wave, "o3", (nt-1)*dt);
        }else{
            sf_putfloat(wave, "d3", snap*dt);
            sf_putfloat(wave, "o3", 0.0f);
        }
    }
        
    /* allocate arrays */
    vv    = sf_floatalloc2(nz, nx);
    dd    = sf_floatalloc2(nt, nx);
    vvpad = sf_floatalloc2(nzpad, nxpad);
    u0    = sf_floatalloc2(nzpad, nxpad);
    u1    = sf_floatalloc2(nzpad, nxpad);
    u2    = sf_floatalloc2(nzpad, nxpad);
    if (pas) {
        mm = NULL;
        ww = sf_floatalloc3(nz, nx, nt);
    } else {
        mm = sf_floatalloc2(nz, nx);
        ww = NULL;
    }
    
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

    memset(u0[0], 0.0f, nzxpad*sizeof(float));
    memset(u1[0], 0.0f, nzxpad*sizeof(float));
    memset(u2[0], 0.0f, nzxpad*sizeof(float));

    /* absorbing boundary condition */
    if (abc) {
        if (verb) sf_warning("absorbing boundary condition");
        abc_init(nzpad,nxpad,nzpad,nxpad,nb,nb,nb,nb,cb,cb,cb,cb);
    }

    if(adj){ /* migration */
        
        /* read data */
        sf_floatread(dd[0], nt*nx, in);
        
        for (it=nt-1; it>-1; it--){
            if (verb) sf_warning("Migration: %d/%d;", it, 0);
            
            /* time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<nxpad-NOP; ix++){
                for (iz=NOP; iz<nzpad-NOP; iz++){
                    u2[ix][iz] = LapT(u1,ix,iz,idx2,idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
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
            for (ix=nb; ix<nb+nx; ix++)
                u1[ix][depth] += dd[ix-nb][it];
            
            if (pas) {
                /* image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
                for     (ix=0; ix<nx; ix++)
                    for (iz=0; iz<nz; iz++)
                        ww[it][ix][iz] = u1[ix+nb][iz+nb];
            }
            if (snap && it%snap==0) sf_floatwrite(u1[0], nzxpad, wave);
        }
        if (verb) sf_warning(".");

        if (!pas) {
            /* output image */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<nx; ix++)
                for (iz=0; iz<nz; iz++)
                    mm[ix][iz] = u1[ix+nb][iz+nb];
            sf_floatwrite(mm[0], nz*nx, out);
        } else {
            /* output source */
            sf_floatwrite(ww[0][0], nz*nx*nt, out);
        }
    
    }else{/* modeling */
    	
        if (pas) { 
            /* read source */
            sf_floatread(ww[0][0], nz*nx*nt, in);
        } else { 
            /* read image */
            sf_floatread(mm[0], nz*nx, in);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<nx; ix++)
                for (iz=0; iz<nz; iz++)
                    u1[ix+nb][iz+nb] = mm[ix][iz];
        }

    	for (it=0; it<nt; it++){
            if (verb) sf_warning("Modeling: %d/%d;", it, nt-1);
            
            if(snap && it%snap==0) sf_floatwrite(u1[0], nzxpad, wave);
            if (pas){
                /* inject source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
                for     (ix=0; ix<nx; ix++)
                    for (iz=0; iz<nz; iz++)
                        u1[ix+nb][iz+nb] += ww[it][ix][iz];
            }

            /* record data */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix=nb; ix<nb+nx; ix++)
	        dd[ix-nb][it] = u1[ix][depth];
            
            if (abc) abc_apply(u0[0]);
            if (abc) abc_apply(u1[0]);
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
        }
    	if (verb) sf_warning(".");
        
    	/* output data */
    	sf_floatwrite(dd[0], nt*nx, out);
    }

    if(pas) {
        free(**ww); free(*ww); free(ww);
    } else {
        free(*mm); free(mm);
    }

    if (abc) abc_close();
    free(*vvpad); free(vvpad); free(*vv); free(vv);
    free(*dd); free(dd); 
    free(*u0); free(u0); free(*u1); free(u1); free(*u2); free(u2);
    exit (0);
}

