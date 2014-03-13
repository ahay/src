/* 2-D prestack reverse time migration and its adjoint for single shot*/
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

static bool verb, snap;
static int nz, nx, nt, nr, is, ns, nw;
static int jt, padx, padz, padnx, padnz;
static int dr, ds, r0, s0, rz, sz, sx;
static float c0, c11, c12, c21, c22;
static float **padvv, *ww;
static sf_file snapshot;

void laplacian(bool adj, float **u0, float **u1, float **u2)
{
    int ix, iz;
    
    if(adj){
#ifdef _OPENMP
#pragma omp parallel for  \
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
#pragma omp parallel for  \
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

void prertm2_oper(bool adj, float **dd, float **mm)
{
    int ix, iz, it, ir, rx, nnt;
    float **u0, **u1, **u2, **temp, **sou2, ***wave;
    
    if(adj)
        memset(mm[0], 0, nz*nx*sizeof(float));
    else
        memset(dd[0], 0, nt*nr*sizeof(float));
    
    nnt=1+(nt-1)/nw;
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    sou2=sf_floatalloc2(nz, nx);
    wave=sf_floatalloc3(nz, nx, nnt);
    
    if(adj){/* migration */
        
        if(verb) sf_warning("ishot/nshot:  %d/%d", is, ns);
            
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
        memset(sou2[0], 0, nz*nx*sizeof(float));
        memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
        
        for(it=0; it<nt; it++){
            sf_warning("RTM_Source: it=%d;",it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
                
            u1[sx][sz]+=ww[it];
            
            if(it%nw==0){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                    wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                }
            }
            }
            if(snap && it%jt==0)
                sf_floatwrite(u1[0], padnx*padnz, snapshot);
        }// end of it
            
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
        
        for(it=nt-1; it>=0; it--){
            sf_warning("RTM_Receiver: it=%d;",it);
            laplacian(false, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
                
            for(ir=0; ir<nr; ir++){
                rx=ir*dr+r0;
                u1[rx][rz]+=dd[ir][it];
            }
            if(it%nw==0){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    mm[ix][iz]+=wave[it/nw][ix][iz]*u1[ix+padx][iz+padz];
                }
            }
            }
        }// end of it
        
        //normalization
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                mm[ix][iz]=mm[ix][iz]/(sou2[ix][iz]+FLT_EPSILON);
        
    }else{/* modeling */
        
        if(verb) sf_warning("ishot/nshot:  %d/%d", is, ns);
        
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
        memset(sou2[0], 0, nz*nx*sizeof(float));
        memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
        
        for(it=0; it<nt; it++){
            sf_warning("Modeling_Source: it=%d;",it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
                
            u1[sx][sz]+=ww[it];
            
            if(it%nw==0){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                    wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                }
            }
            }
            if(snap && it%jt==0)
                sf_floatwrite(u1[0], padnx*padnz, snapshot);
        }// end of it
            
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
        
        //normalization
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                mm[ix][iz]=mm[ix][iz]/(sou2[ix][iz]+FLT_EPSILON);
        
        for(it=0; it<nt; it++){
            sf_warning("Modeling_Receiver: it=%d;",it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
            
            if(it%nw==0){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    u1[ix+padx][iz+padz]+=wave[it/nw][ix][iz]*mm[ix][iz];
                }
            }
            }
            for(ir=0; ir<nr; ir++){
                rx=ir*dr+r0;
                dd[ir][it]+=u1[rx][rz];
            }
        } //end of it
    }// end of if
}

int main(int argc, char* argv[])
{
    bool adj;
    int ix, iz;
    float dz, dx, dt, z0, x0, t0;
    double dt2, idx2, idz2;
 
    float **dd, **mm, **vv;
    sf_file in, out, vel, wavelet;
    
    sf_init(argc, argv);
    omp_init();
    
    if(!sf_getbool("adj", &adj)) adj=true;
    if(!sf_getbool("verb", &verb)) verb=true;
    if(!sf_getbool("snap", &snap)) snap=false;
    
    in=sf_input("in");
    out=sf_output("out");
    vel=sf_input("velocity");
    wavelet=sf_input("wavelet");
    
    if(!sf_histint(vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "o1", &z0)) sf_error("No o1= in velocity");
    
    if(!sf_histint(vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");
    if(!sf_histfloat(vel, "o2", &x0)) sf_error("No o2= in velocity");
    
    if(!sf_histint(wavelet, "n1", &nt)) sf_error("No n1= in wavelet");
    if(!sf_histfloat(wavelet, "d1", &dt)) sf_error("No d1= in wavelet");
    if(!sf_histfloat(wavelet, "o1", &t0)) sf_error("No o1= in wavelet");
    
    if(adj){/* migration */
        if(!sf_histint(in, "n2", &nr)) sf_error("No n2= in input");
        if(!sf_histint(in, "d2", &dr)) sf_error("No d2= in input");
        if(!sf_histint(in, "o2", &r0)) sf_error("No o2= in input");
        
        sf_putint(out, "n1", nz);
        sf_putfloat(out, "d1", dz);
        sf_putfloat(out, "o1", z0);
        sf_putint(out, "n2", nx);
        sf_putfloat(out, "d2", dx);
        sf_putfloat(out, "o2", x0);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "km");
        sf_putstring(out, "label2", "Lateral");
        sf_putstring(out, "unit2", "km");
    }else{/* modeling */
        if(!sf_getint("nr", &nr)) nr=nx;
        if(!sf_getint("dr", &dr)) dr=1;
        if(!sf_getint("r0", &r0)) r0=0;
        
        sf_putint(out, "n1", nt);
        sf_putfloat(out, "d1", dt);
        sf_putfloat(out, "o1", t0);
        sf_putint(out, "n2", nr);
        sf_putint(out, "d2", dr);
        sf_putint(out, "o2", r0);
        sf_putstring(out, "label1", "Time");
        sf_putstring(out, "unit1", "s");
        sf_putstring(out, "label2", "Reciver");
        sf_putstring(out, "unit2", "");
    }
    
    if(!sf_getint("is", &is)) sf_error("Need is=");
    /* shot index, starting from 0 */
    if(!sf_getint("ns", &ns)) sf_error("Need ns=");
    /* shot number */
    if(!sf_getint("ds", &ds)) sf_error("Need ds=");
    if(!sf_getint("s0", &s0)) s0=0;
    
    if(!sf_getint("nw", &nw)) nw=1;
    
    if(!sf_getint("rz", &rz)) rz=0;
    if(!sf_getint("sz", &sz)) sz=0;
    if(!sf_getint("jt", &jt)) jt=100;
    if(!sf_getint("padz", &padz)) padz=nz;
    if(!sf_getint("padx", &padx)) padx=nz;
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    rz=padz+rz;
    sz=padz+sz;
    s0=s0+padx;
    sx=s0+is*ds;
    r0=r0+padx;
    
    dd=sf_floatalloc2(nt, nr);
    mm=sf_floatalloc2(nz, nx);
    vv=sf_floatalloc2(nz, nx);
    padvv=sf_floatalloc2(padnz, padnx);
    ww=sf_floatalloc(nt);
    
    sf_floatread(ww, nt, wavelet);
    sf_floatread(vv[0], nz*nx, vel);
    
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
    
    if(snap){
        snapshot=sf_output("snapshot");
        
        sf_putint(snapshot, "n1", padnz);
        sf_putint(snapshot, "d1", 1);
        sf_putint(snapshot, "o1", -padz);
        sf_putint(snapshot, "n2", padnx);
        sf_putint(snapshot, "d2", 1);
        sf_putint(snapshot, "o2", -padx);
        sf_putint(snapshot, "n3", 1+(nt-1)/jt);
        sf_putfloat(snapshot, "d3", jt*dt);
        sf_putfloat(snapshot, "o3", t0);
        sf_putstring(snapshot, "label1", "Depth");
        sf_putstring(snapshot, "unit1", "");
        sf_putstring(snapshot, "label2", "Lateral");
        sf_putstring(snapshot, "unit2", "");
        sf_putstring(snapshot, "label3", "Time");
        sf_putstring(snapshot, "unit3", "s");
    }
    
    sf_warning("nx=%d padnx=%d nz=%d padnz=%d nt=%d nr=%d sx=%d", nx, padnx, nz, padnz, nt, nr, ns);
    sf_warning("dx=%.3f dz=%.3f dt=%.3f dr=%d ds=%d", dx, dz, dt, dr, ds);
    sf_warning("x0=%.3f z0=%.3f t0=%.3f r0=%d s0=%d", x0, z0, t0, r0, s0);
    
    idz2=1./(dz*dz);
    idx2=1./(dx*dx);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2*(c11+c12+c21+c22);
    
    if(adj){/* migration */
        sf_floatread(dd[0], nt*nr, in);
        prertm2_oper(adj, dd, mm);
        sf_floatwrite(mm[0], nz*nx, out);
    }else{/* modeling */
        sf_floatread(mm[0], nz*nx, in);
        prertm2_oper(adj, dd, mm);
        sf_floatwrite(dd[0], nt*nr, out);
    }
    
    exit(0);
}
