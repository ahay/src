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
static int nz, nx, nt, nr, ns;
static int jt, padx, padz, padnx, padnz;
static int dr_v, ds_v, r0_v, s0_v, zr_v, zs_v;
static float c0, c11, c12, c21, c22;
static double idx2, idz2;
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

void prertm2_oper(bool adj, int is, float **localdd, float **localmm)
{
    int ix, iz, it, ir, sx, rx;
    float **u0, **u1, **u2, **temp, **sou2, ***wave;
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    sou2=sf_floatalloc2(nz, nx);
    wave=sf_floatalloc3(nz, nx, nt);
    
    memset(u0[0], 0, padnz*padnx*sizeof(float));
    memset(u1[0], 0, padnz*padnx*sizeof(float));
    memset(u2[0], 0, padnz*padnx*sizeof(float));
    memset(sou2[0], 0, nz*nx*sizeof(float));
    memset(wave[0][0], 0, nz*nx*nt*sizeof(float));
    
    sx=is*ds_v+s0_v;
    
    if(adj){/* migration */
        
        for(it=0; it<nt; it++){
            if(verb) sf_warning("is=%d; RTM_Source: it=%d;", is+1, it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
            
            u1[sx][zs_v]+=ww[it];
            
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                    wave[it][ix][iz]=u1[ix+padx][iz+padz];
                }
            }
                
            if(snap && is==0 && it%jt==0)
                sf_floatwrite(u1[0], padnx*padnz, snapshot);
        }// end of it
            
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
            
        for(it=nt-1; it>=0; it--){
            if(verb) sf_warning("is=%d; RTM_Receiver: it=%d;", is+1, it);
            laplacian(false, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
                
            for(ir=0; ir<nr; ir++){
                rx=ir*dr_v+r0_v;
                u1[rx][zr_v]+=localdd[ir][it];
            }

            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    localmm[ix][iz]+=wave[it][ix][iz]*u1[ix+padx][iz+padz];
                }
            }
        }//end of it
        
        /* normalization 
        for(ix=0; ix<nx; ix++){
            for(iz=0; iz<nz; iz++){
                localmm[ix][iz]/=(sou2[ix][iz]+FLT_EPSILON);
            }
        } */
        
    }else{/* modeling */
            
        for(it=0; it<nt; it++){
            if(verb) sf_warning("is=%d; Modeling_Source: it=%d;",is+1, it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
            
            u1[sx][zs_v]+=ww[it];
            
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                    wave[it][ix][iz]=u1[ix+padx][iz+padz];
                }
            }
                
            if(snap && is==0 && it%jt==0)
                sf_floatwrite(u1[0], padnx*padnz, snapshot);
        }// end of it
            
        memset(u0[0], 0, padnz*padnx*sizeof(float));
        memset(u1[0], 0, padnz*padnx*sizeof(float));
        memset(u2[0], 0, padnz*padnx*sizeof(float));
        
        /* normalization 
        for(ix=0; ix<nx; ix++){
            for(iz=0; iz<nz; iz++){
                localmm[ix][iz]/=(sou2[ix][iz]+FLT_EPSILON);
         }
         } */
            
        for(it=0; it<nt; it++){
            if(verb) sf_warning("is=%d; Modeling_Receiver: it=%d;",is+1, it);
            laplacian(true, u0, u1, u2);
                
            temp=u0;
            u0=u1;
            u1=u2;
            u2=temp;
            
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    u1[ix+padx][iz+padz]+=wave[it][ix][iz]*localmm[ix][iz];
                }
            }
                
            for(ir=0; ir<nr; ir++){
                rx=ir*dr_v+r0_v;
                localdd[ir][it]+=u1[rx][zr_v];
            }
        } //end of it
    }// end of if
}

int main(int argc, char* argv[])
{
    bool adj;
    int ix, iz, is, ir, it, iturn;
    int nspand, cpuid, numprocs;
    float dz, dx, dt, dr, ds;
    float z0, x0, t0, r0, s0;
    float zr, zs, padx0, padz0;
    float dt2;
    float *sendbuf, *recvbuf;
    double tstart, tend, duration, *sendbuf2;
 
    float ***dd, **localdd, **mm, **localmm, **vv;
    sf_file in, out, vel, wavelet;
    
    sf_init(argc, argv);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    tstart=MPI_Wtime();
    
    if(cpuid==0) sf_warning("numprocs=%d", numprocs);
    
    if(!sf_getbool("adj", &adj)) adj=true;
    if(!sf_getbool("verb", &verb)) verb=false;
    if(!sf_getbool("snap", &snap)) snap=false;
    
    in=sf_input("--input");
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
    }else{/* modeling */
        if(!sf_getint("nr", &nr)) sf_error("No nr=");
        if(!sf_getfloat("dr", &dr)) sf_error("No dr=");
        if(!sf_getfloat("r0", &r0)) sf_error("No r0=");
        
        if(!sf_getint("ns", &ns)) sf_error("No ns=");
        if(!sf_getfloat("ds", &ds)) sf_error("No ds=");
        if(!sf_getfloat("s0", &s0)) sf_error("No s0=");
    }
    
    if(ns%numprocs==0)
        nspand=ns;
    else
        nspand=(ns/numprocs+1)*numprocs;
    
    if(cpuid==0){
        out=sf_output("--output");
        if(adj){
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
        }else{
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
    }
    
    if(!sf_getfloat("zr", &zr)) zr=0.0;
    if(!sf_getfloat("zs", &zs)) zs=0.0;
    if(!sf_getint("jt", &jt)) jt=100;
    if(!sf_getint("padz", &padz)) padz=nz;
    if(!sf_getint("padx", &padx)) padx=nz;
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    padx0=x0-dx*padx;
    padz0=z0-dz*padz;
    
    if(cpuid==0) dd=sf_floatalloc3(nt, nr, nspand);
    mm=sf_floatalloc2(nz, nx);
    localdd=sf_floatalloc2(nt, nr);
    localmm=sf_floatalloc2(nz, nx);
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
    
    if(snap && cpuid==0){
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
    
    if(cpuid==0){
    sf_warning("nx=%d padnx=%d nz=%d padnz=%d nt=%d nr=%d ns=%d", nx, padnx, nz, padnz, nt, nr, ns);
    sf_warning("dx=%.3f dz=%.3f dt=%.3f dr=%.3f ds=%.3f", dx, dz, dt, dr, ds);
    sf_warning("x0=%.3f z0=%.3f t0=%.3f r0=%.3f s0=%.3f", x0, z0, t0, r0, s0);
    sf_warning("dr_v=%d r0_v=%d zr_v=%d ds_v=%d s0_v=%d zs_v=%d", dr_v, r0_v, zr_v, ds_v, s0_v, zs_v);
    }
    
    idz2=1./(dz*dz);
    idx2=1./(dx*dx);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2*(c11+c12+c21+c22);
    
    if(adj){/* migration */
        memset(mm[0], 0., nz*nx*sizeof(float));
        
        if(cpuid==0){
            sf_floatread(dd[0][0], nt*nr*ns, in);
            for(is=ns; is<nspand; is++)
                for(ir=0; ir<nr; ir++)
                    for(it=0; it<nt; it++)
                        dd[is][ir][it]=0.0;
        }
        
        for(iturn=0; iturn*numprocs<nspand; iturn++){
            is=iturn*numprocs+cpuid;
            
            sf_warning("********* Reverse Time Migration Process ***********\n");
            sf_warning("Processor ID: %d      Current shot number: %d\n", cpuid, is+1);
            sf_warning("****************************************************\n");
            
            if(cpuid==0){
                sendbuf=dd[iturn*numprocs][0];
                recvbuf=localdd[0];
            }else{
                sendbuf=NULL;
                recvbuf=localdd[0];
            }
            MPI_Scatter(sendbuf, nr*nt, MPI_FLOAT, recvbuf, nr*nt, MPI_FLOAT, 0, MPI_COMM_WORLD);
            
            memset(localmm[0], 0., nz*nx*sizeof(float));
            if(is<ns) prertm2_oper(adj, is, localdd, localmm);
            
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    mm[ix][iz]+=localmm[ix][iz];
        }
        
        if(cpuid==0){
            sendbuf=MPI_IN_PLACE;
            recvbuf=mm[0];
        }else{
            sendbuf=mm[0];
            recvbuf=NULL;
        }
        MPI_Reduce(sendbuf, recvbuf, nz*nx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(cpuid==0) sf_floatwrite(mm[0], nz*nx, out);
        
    }else{ /* modeling */
        if(cpuid==0) sf_floatread(mm[0], nz*nx, in);
        
        MPI_Bcast(mm[0], nz*nx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        for(iturn=0; iturn*numprocs<nspand; iturn++){
            is=iturn*numprocs+cpuid;
            
            sf_warning("**** Forward Time Modeling Process ******\n");
            sf_warning("Processor ID: %d  Current shot number: %d\n", cpuid, is+1);
            sf_warning("*****************************************\n");
            
            memset(localdd[0], 0., nr*nt*sizeof(float));
            if(is<ns) prertm2_oper(adj, is, localdd, mm);
            
            if(cpuid==0){
                sendbuf=localdd[0];
                recvbuf=dd[iturn*numprocs][0];
            }else{
                sendbuf=localdd[0];
                recvbuf=NULL;
            }
            MPI_Gather(sendbuf, nt*nr, MPI_FLOAT, recvbuf, nt*nr, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        
        if(cpuid==0) sf_floatwrite(dd[0][0], ns*nr*nt, out);
    }
    
    sf_fileclose(in); sf_fileclose(vel); sf_fileclose(wavelet);
    if(cpuid==0){
        if(snap) sf_fileclose(snapshot);
        sf_fileclose(out);
        free(**dd); free(*dd); free(dd);
    }
    free(*padvv); free(padvv);
    free(ww);
    free(*vv); free(vv);
    free(*mm); free(mm);
    free(*localdd); free(localdd);
    free(*localmm); free(localmm);
    
    tend=MPI_Wtime();
    duration=tend-tstart;
    if(cpuid==0) sendbuf2=MPI_IN_PLACE;
    else sendbuf2=&duration;
    MPI_Reduce(sendbuf2, &duration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(cpuid==0) sf_warning("The running time: %e\n", duration);
    
    MPI_Finalize();
    
    exit(0);
}
