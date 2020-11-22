/* 2-D exploding-reflector RTM and its adjoint */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

static int padnx, padnz;
static float c0, c11, c12, c21, c22;
static float **padvv;

void laplacian(bool adj, float **u0, float **u1, float **u2)
{
    int ix, iz;
    
    if(adj){
#ifdef _OPENMP
#pragma omp parallel for  default(none) \
private(ix, iz)       \
shared(padnx, padnz, u0, u1, u2, padvv, c0, c11, c12, c21, c22)
#endif
        for(ix=2; ix<padnx-2; ix++){
            for(iz=2; iz<padnz-2; iz++){
                u2[ix][iz]=
                (c11*(u1[ix][iz-1]+u1[ix][iz+1])+
                 c12*(u1[ix][iz-2]+u1[ix][iz+2])+
                 c0*u1[ix][iz]+
                 c21*(u1[ix-1][iz]+u1[ix+1][iz])+
                 c22*(u1[ix-2][iz]+u1[ix+2][iz]))
                *padvv[ix][iz]+2*u1[ix][iz]-u0[ix][iz];
            }
        }
    }else{
#ifdef _OPENMP
#pragma omp parallel for  default(none) \
private(ix, iz)       \
shared(padnx, padnz, u0, u1, u2, padvv, c0, c11, c12, c21, c22)
#endif
        for(ix=2; ix<padnx-2; ix++){
            for(iz=2; iz<padnz-2; iz++){
                u2[ix][iz]=
                (c11*(u1[ix][iz-1]*padvv[ix][iz-1]+u1[ix][iz+1]*padvv[ix][iz+1])+
                 c12*(u1[ix][iz-2]*padvv[ix][iz-2]+u1[ix][iz+2]*padvv[ix][iz+2])+
                 c0*u1[ix][iz]*padvv[ix][iz]+
                 c21*(u1[ix-1][iz]*padvv[ix-1][iz]+u1[ix+1][iz]*padvv[ix+1][iz])+
                 c22*(u1[ix-2][iz]*padvv[ix-2][iz]+u1[ix+2][iz]*padvv[ix+2][iz]))
                +2*u1[ix][iz]-u0[ix][iz];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    bool adj, snap; /* adjoint flag */
    int ix, iz, it; /* index variables */
    int nt, nx, nz, n0, jt, n12, padx, padz, n2;
    float dt, dx, dz, dt2, idz2, idx2;
    
    float **dd, **mm, **vv;
    float **u0, **u1, **u2, **tmp; /* temporary arrays */

    sf_file in, out, vel, wave=NULL; /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);
    
    /* initialize OpenMP support */
#ifdef _OPENMP
    omp_init();
#endif

    if(!sf_getbool("adj", &adj)) adj=true;
    /* adjoint flag, 0: modeling, 1: migration */
    if(!sf_getbool("snap", &snap)) snap=false;
    /* wavefield snapshot flag */
    if(!sf_getint("n0", &n0)) n0=0;
    /* surface */
    if(!sf_getint("jt", &jt)) jt=50;
    /* time interval of wavefield snapshot */
	
    /* setup I/O files */
    in=sf_input("in");
    out=sf_output("out");
    vel=sf_input("velocity");
    /* velocity model */
    
    /* Dimensions */
    if(!sf_histint(vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histint(vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");

    if(adj){ /* migration */
        if(!sf_histint(in, "n1", &nt)) sf_error("No n1= in data");
        if(!sf_histfloat(in, "d1", &dt)) sf_error("No d1= in data");
        if(!sf_histint(in, "n2", &n2) || n2!=nx) sf_error("Need n2=%d in data", nx);
        
        sf_putint(out, "n1", nz);
        sf_putfloat(out, "d1", dz);
        sf_putfloat(out, "o1", 0.0);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "km");
        sf_putstring(out, "label2", "Lateral");
        sf_putstring(out, "unit2", "km");
    }else{ /* modeling */
        if(!sf_getint("nt", &nt)) sf_error("Need nt=");
        if(!sf_getfloat("dt", &dt)) sf_error("Need dt=");
        
        sf_putint(out, "n1", nt);
        sf_putfloat(out, "d1", dt);
        sf_putfloat(out, "o1", 0.0);
        sf_putstring(out, "label1", "Time");
        sf_putstring(out, "unit1", "s");
        sf_putstring(out, "label2", "Lateral");
        sf_putstring(out, "unit2", "km");
    }
    
    /* lengths of padding boundary */
    if(!sf_getint("padx", &padx)) padx=nz/2;
    if(!sf_getint("padz", &padz)) padz=nz/2;
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    n0=n0+padz;
    n12=padnz*padnx;

    /* set Laplacian coefficients */
    idz2=1.0/(dz*dz);
    idx2=1.0/(dx*dx);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2.0*(c11+c12+c21+c22);
    
    /* wavefield snapshot */
    if(snap){
        wave=sf_output("wave");
        
        sf_putint(wave, "n1", padnz);
        sf_putint(wave, "d1", 1);
        sf_putint(wave, "o1", -padz);
            
        sf_putint(wave, "n2", padnx);
        sf_putint(wave, "d2", 1);
        sf_putint(wave, "o2", -padx);
        
        sf_putint(wave, "n3", 1+(nt-1)/jt);
        if(adj){
            sf_putfloat(wave, "d3", -jt*dt);
            sf_putfloat(wave, "o3", (nt-1)*dt);
        }else{
            sf_putfloat(wave, "d3", jt*dt);
            sf_putfloat(wave, "o3", 0.0);
        }
    }
        
    /* allocate arrays */
    vv=sf_floatalloc2(nz, nx);
    dd=sf_floatalloc2(nt, nx);
    mm=sf_floatalloc2(nz, nx);
    padvv=sf_floatalloc2(padnz, padnx);
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    
    /* read velocity */
    sf_floatread(vv[0], nz*nx, vel);
    
    /* pad boundary */
    dt2=dt*dt;
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix+padx][iz+padz]=vv[ix][iz]*vv[ix][iz]*dt2;
    for(iz=0; iz<padz; iz++){
        for(ix=padx; ix<nx+padx; ix++){
            padvv[ix][iz]=padvv[ix][padz];
            padvv[ix][iz+nz+padz]=padvv[ix][nz+padz-1];
        }
    }
    for(ix=0; ix<padx; ix++){
        for(iz=0; iz<padnz; iz++){
            padvv[ix][iz]=padvv[padx][iz];
            padvv[ix+nx+padx][iz]=padvv[nx+padx-1][iz];
        }
    }
    
    memset(u0[0], 0.0, n12*sizeof(float));
    memset(u1[0], 0.0, n12*sizeof(float));
    memset(u2[0], 0.0, n12*sizeof(float));

    if(adj){ /* migration */
        
        /* read data */
        sf_floatread(dd[0], nt*nx, in);
        
        for(it=nt-1; it>=0; it--){
            sf_warning("Migration: %d;", it);
            
            laplacian(adj, u0, u1, u2);
            tmp=u0; u0=u1; u1=u2; u2=tmp;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix) shared(padx, nx, dd, u1, it, n0)
#endif
            for(ix=padx; ix<padx+nx; ix++)
            /* inject data */
            u1[ix][n0]+=dd[ix-padx][it];
            
            if(snap && it%jt==0) sf_floatwrite(u1[0], n12, wave);
        }
        sf_warning(".");
        
        /* output image */
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                mm[ix][iz]=u1[ix+padx][iz+padz];
        sf_floatwrite(mm[0], nz*nx, out);
    
    }else{/* modeling */
    	
        /* read reflector */
    	sf_floatread(mm[0], nz*nx, in);
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                u1[ix+padx][iz+padz]=mm[ix][iz];
    	    
    	for(it=0; it<nt; it++){
            sf_warning("Modeling: %d;", it);
            
            if(snap && it%jt==0) sf_floatwrite(u1[0], n12, wave);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix) shared(padx, nx, dd, u1, it, n0)
#endif
            for(ix=padx; ix<padx+nx; ix++)
            /* record data */
	        dd[ix-padx][it]=u1[ix][n0];
            
            laplacian(adj, u0, u1, u2);
            tmp=u0; u0=u1; u1=u2; u2=tmp;
        }
    	sf_warning(".");
        
    	/* output data */
    	sf_floatwrite(dd[0], nt*nx, out);
    }
    
    free(*padvv); free(padvv); free(*vv); free(vv);
    free(*dd); free(dd); free(*mm); free(mm);
    free(*u0); free(u0); free(*u1); free(u1); free(*u2); free(u2);
    exit (0);
}
