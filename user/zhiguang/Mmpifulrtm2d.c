/* 2-D prestack reverse time migration and its adjoint with MPI for full coverage*/
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
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static bool verb, snap;
static int nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay;
static int jt, padx, padz, padnx, padnz;
static int dr_v, ds_v, r0_v, s0_v, zr_v, zs_v;
static float c0, c11, c12, c21, c22;
static double idx2, idz2;
static float **padvv, *ww, tdelay;
static sf_file in, out, snapshot;
static int cpuid, numprocs;

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

void prertm2_oper(bool adj, float **mm)
{
    int ix, iz, is, it, ir, sx, rx, nnt;
    float **u0, **u1, **u2, **temp, **sou2, **perm, ***wave, **localsum, **dd;
    FILE *swap;
    
    swap=fopen("temswap.bin","wb+");
    if(adj) memset(mm[0], 0, nz*nx*sizeof(float));
    nnt=1+(nt-1)/nw;
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    sou2=sf_floatalloc2(nz, nx);
    perm=sf_floatalloc2(nz, nx);
    localsum=sf_floatalloc2(nz, nx);
    wave=sf_floatalloc3(nz, nx, nnt);
    dd=sf_floatalloc2(nt, nr);
    memset(localsum[0], 0, nz*nx*sizeof(float));
    
    if(adj){/* migration */
        
        for(is=cpuid; is<ns; is+=numprocs){
            if(verb) sf_warning("ishot @@@ nshot:%d  %d cpuid=%d migration", is+1, ns, cpuid);
            sx=is*ds_v+s0_v;
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(perm[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
            
            sf_seek(in, is*nr*nt*sizeof(float), SEEK_SET);
            sf_floatread(dd[0], nr*nt, in);
            
            for(it=0; it<nt; it++){
                //sf_warning("RTM_Source: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ix=0; ix<nsource; ix++){
                if(it>=ix*ndelay)
                u1[sx+ix*dsource][zs_v]+=ww[it-ix*ndelay];
                }
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                        wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                    }
                }
                }
                
                if(snap && is==ns/2 && it%jt==0)
                    sf_floatwrite(u1[0], padnx*padnz, snapshot);
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            
            for(it=nt-1; it>=0; it--){
                //sf_warning("RTM_Receiver: it=%d;",it);
                laplacian(false, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ir=0; ir<nr; ir++){
                    rx=ir*dr_v+r0_v;
                    u1[rx][zr_v]+=dd[ir][it];
                }
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        perm[ix][iz]+=wave[it/nw][ix][iz]*u1[ix+padx][iz+padz];
                    }
                }
                }
            }// end of it
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    localsum[ix][iz]+=perm[ix][iz]/(sou2[ix][iz]+FLT_EPSILON);
        }// end of shot
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(localsum[0], mm[0], nz*nx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        
    }else{/* modeling */
        
        for(is=cpuid; is<ns; is+=numprocs){
            if(verb) sf_warning("ishot @@@ nshot:%d  %d cpuid=%d modeling", is+1, ns, cpuid);
            sx=is*ds_v+s0_v;
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(perm[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
            memset(dd[0], 0, nt*nr*sizeof(float));
            
            for(it=0; it<nt; it++){
                //sf_warning("RTM_Source: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ix=0; ix<nsource; ix++){
                if(it>=ix*ndelay)
                u1[sx+ix*dsource][zs_v]+=ww[it-ix*ndelay];
                }
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                        wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                    }
                }
                }
                
                if(snap && is==ns/2 && it%jt==0)
                    sf_floatwrite(u1[0], padnx*padnz, snapshot);
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    perm[ix][iz]+=mm[ix][iz]/(sou2[ix][iz]+FLT_EPSILON);
            
            for(it=0; it<nt; it++){
                //sf_warning("Modeling_Receiver: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        u1[ix+padx][iz+padz]+=wave[it/nw][ix][iz]*perm[ix][iz];
                    }
                }
                }
                
                for(ir=0; ir<nr; ir++){
                    rx=ir*dr_v+r0_v;
                    dd[ir][it]+=u1[rx][zr_v];
                }
            } //end of it
            
            fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
            fwrite(dd[0], sizeof(float)*nt, nr, swap);
        }// end of shot
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(cpuid==0){
            for(is=0; is<ns; is++){
                fseeko(swap, is*nr*nt*sizeof(float), SEEK_SET);
                fread(dd[0], sizeof(float)*nt, nr, swap);
                sf_floatwrite(dd[0], nr*nt, out);
            }
            fclose(swap);
            remove("temswap.bin");
        }
    MPI_Barrier(MPI_COMM_WORLD);
    }// end of if
}

int main(int argc, char* argv[])
{
    bool adj;
    int ix, iz;
    float dz, dx, dt, dr, ds;
    float z0, x0, t0, r0, s0;
    float zr, zs, padx0, padz0;
    float dt2;
 
    float ***dd, **mm, **vv;
    sf_file vel, wavelet;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    sf_init(argc, argv);
//    omp_init();
    
    if(cpuid==0) sf_warning("numprocs=%d", numprocs);
    
    if(!sf_getbool("adj", &adj)) adj=true;
    if(!sf_getbool("verb", &verb)) verb=true;
    if(!sf_getbool("snap", &snap)) snap=false;
    
    in=sf_input("input");
    out=sf_output("output");
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
        if(!sf_histfloat(in, "d2", &dr)) sf_error("No d2= in input");
        if(!sf_histfloat(in, "o2", &r0)) sf_error("No o2= in input");
        
        if(!sf_histint(in, "n3", &ns)) sf_error("No n3= in input");
        if(!sf_histfloat(in, "d3", &ds)) sf_error("No d3= in input");
        if(!sf_histfloat(in, "o3", &s0)) sf_error("No o3= in input");
        
        sf_putint(out, "n1", nz);
        sf_putfloat(out, "d1", dz);
        sf_putfloat(out, "o1", z0);
        sf_putint(out, "n2", nx);
        sf_putfloat(out, "d2", dx);
        sf_putfloat(out, "o2", x0);
        sf_putint(out, "n3", 1);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "km");
        sf_putstring(out, "label2", "Lateral");
        sf_putstring(out, "unit2", "km");
    }else{/* modeling */
        if(!sf_getint("nr", &nr)) sf_error("No nr=");
        if(!sf_getfloat("dr", &dr)) sf_error("No dr=");
        if(!sf_getfloat("r0", &r0)) sf_error("No r0=");
        
        if(!sf_getint("ns", &ns)) sf_error("No ns=");
        if(!sf_getfloat("ds", &ds)) sf_error("No ds=");
        if(!sf_getfloat("s0", &s0)) sf_error("No s0=");
        
        sf_putint(out, "n1", nt);
        sf_putfloat(out, "d1", dt);
        sf_putfloat(out, "o1", t0);
        sf_putint(out, "n2", nr);
        sf_putfloat(out, "d2", dr);
        sf_putfloat(out, "o2", r0);
        sf_putint(out, "n3", ns);
        sf_putfloat(out, "d3", ds);
        sf_putfloat(out, "o3", s0);
        sf_putstring(out, "label1", "Time");
        sf_putstring(out, "unit1", "s");
        sf_putstring(out, "label2", "Lateral");
        sf_putstring(out, "unit2", "km");
        sf_putstring(out, "label3", "Shot");
        sf_putstring(out, "unit3", "km");
    }
    
    if(!sf_getfloat("zr", &zr)) zr=0.0;
    if(!sf_getfloat("zs", &zs)) zs=0.0;
    if(!sf_getint("jt", &jt)) jt=100;
    if(!sf_getint("padz", &padz)) padz=nz;
    if(!sf_getint("padx", &padx)) padx=nz;
    if(!sf_getint("nw", &nw)) nw=1;
    if(!sf_getint("nsource", &nsource)) nsource=1;
    if(!sf_getint("dsource", &dsource)) dsource=0;
    if(!sf_getfloat("tdelay", &tdelay)) tdelay=0;
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    padx0=x0-dx*padx;
    padz0=z0-dz*padz;
    
    dd=sf_floatalloc3(nt, nr, ns);
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
        sf_putfloat(snapshot, "d1", dz);
        sf_putfloat(snapshot, "o1", padz0);
        sf_putint(snapshot, "n2", padnx);
        sf_putfloat(snapshot, "d2", dx);
        sf_putfloat(snapshot, "o2", padx0);
        sf_putint(snapshot, "n3", 1+(nt-1)/jt);
        sf_putfloat(snapshot, "d3", jt*dt);
        sf_putfloat(snapshot, "o3", t0);
        sf_putstring(snapshot, "label1", "Depth");
        sf_putstring(snapshot, "unit1", "km");
        sf_putstring(snapshot, "label2", "Lateral");
        sf_putstring(snapshot, "unit2", "km");
        sf_putstring(snapshot, "label3", "Time");
        sf_putstring(snapshot, "unit3", "s");
    }
    
    dr_v=(dr/dx)+0.5;
    r0_v=(r0-padx0)/dx+0.5;
    zr_v=(zr-padz0)/dz+0.5;
    
    ds_v=(ds/dx)+0.5;
    s0_v=(s0-padx0)/dx+0.5;
    zs_v=(zs-padz0)/dz+0.5;
    ndelay=tdelay/dt;
    
    if(cpuid==0){
    sf_warning("nx=%d padnx=%d nz=%d padnz=%d nt=%d nr=%d ns=%d", nx, padnx, nz, padnz, nt, nr, ns);
    sf_warning("dx=%.3f dz=%.3f dt=%.3f dr=%.3f ds=%.3f", dx, dz, dt, dr, ds);
    sf_warning("x0=%.3f z0=%.3f t0=%.3f r0=%.3f s0=%.3f", x0, z0, t0, r0, s0);
    sf_warning("dr_v=%d r0_v=%d zr_v=%d ds_v=%d s0_v=%d zs_v=%d", dr_v, r0_v, zr_v, ds_v, s0_v, zs_v);
    sf_warning("nsource=%d dsource=%d tdelay=%.3f ndelay=%d", nsource, dsource, tdelay, ndelay);
    }
    
    idz2=1./(dz*dz);
    idx2=1./(dx*dx);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2*(c11+c12+c21+c22);
    
    if(adj){/* migration */
        prertm2_oper(adj, mm);
        if(cpuid==0) sf_floatwrite(mm[0], nz*nx, out);
    }else{
        sf_floatread(mm[0], nz*nx, in);
        prertm2_oper(adj, mm);
    }
    
    sf_fileclose(in); sf_fileclose(out);
    free(*padvv); free(padvv); free(ww);
    MPI_Finalize();
    
    exit(0);
}
