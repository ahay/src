/* 2-D prestack reverse-time migration and its adjoint */
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
static int nz, nx, nvx, nt, nr, ns;
static int jt, padx, padz, padnx, padnz;
static int dr_v, r0_v, zr_v, ds_v, sx, zs_v;
static float c0, c11, c12, c21, c22;
static float **vv, **padvv, *ww;
static sf_file snapshot, in, out;
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

void get_padvv(int is)
{
    int ix, iz;
    int v1, v2;
    
    v1=is*ds_v;
    v2=v1+nx;
    for(ix=v1; ix<v2; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix-v1+padx][iz+padz]=vv[ix][iz];
    
    for(ix=0; ix<padx; ix++){
        for(iz=padz; iz<nz+padz; iz++){
            padvv[ix][iz]=padvv[padx][iz];
            padvv[nx+padx+ix][iz]=padvv[nx+padx-1][iz];
        }
    }
    
    for(iz=0; iz<padz; iz++){
        for(ix=0; ix<padnx; ix++){
            padvv[ix][iz]=padvv[ix][padz];
            padvv[ix][nz+padz+iz]=padvv[ix][nz+padz-1];
        }
    }
}

void prertm2_oper(bool adj, float **mm)
{
    int ix, iz, is, it, ir, rx, v1, v2;
    float **u0, **u1, **u2, **temp, **permm, **sou2, ***wave, **localsum;
    float **dd;
    FILE *swap;
    
    swap=fopen("temswap.bin","wb");
    
    if(adj) memset(mm[0], 0, nz*nvx*sizeof(float));
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    permm=sf_floatalloc2(nz, nx);
    sou2=sf_floatalloc2(nz, nx);
    wave=sf_floatalloc3(nz, nx, nt);
    dd=sf_floatalloc2(nt, nr);
    localsum=sf_floatalloc2(nz, nvx);
    
    if(adj){/* migration */
        
        for(is=cpuid; is<ns; is+=numprocs){
            if(verb) sf_warning("ishot @@@ nshot: %d %d cupid=%d migration", is+1, ns, cpuid);
            
            get_padvv(is);
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nt*sizeof(float));
            
            sf_seek(in, is*nr*nt*sizeof(float), SEEK_SET);
            sf_floatread(dd[0], nr*nt, in);
            
            for(it=0; it<nt; it++){
                //sf_warning("RTM_Source: it=%d;",it);
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
                if(snap && is< numprocs && cpuid==0 && it%jt==0)
                    sf_floatwrite(u1[0], padnx*padnz, snapshot);
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(permm[0], 0, nz*nx*sizeof(float));
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
                
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        permm[ix][iz]+=wave[it][ix][iz]*u1[ix+padx][iz+padz];
                    }
                }
            }// end of it
            
            v1=is*ds_v;
            v2=v1+nx;
            for(ix=v1; ix<v2; ix++){
                for(iz=0; iz<nz; iz++){
                    localsum[ix][iz]+=permm[ix-v1][iz]/(sou2[ix-v1][iz]+FLT_EPSILON);
                }
            }                        
        }// end of shot
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(localsum[0], mm[0], nz*nvx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        
    }else{/* modeling */
        
        for(is=cpuid; is<ns; is+=numprocs){
            if(verb) sf_warning("ishot @@@ nshot: %d %d cpuid=%d modeling", is+1, ns, cpuid);
            
            get_padvv(is);
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nt*sizeof(float));
            memset(dd[0], 0, nr*nt*sizeof(float));
            for(it=0; it<nt; it++){
                //sf_warning("Modeling_Source: it=%d;",it);
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
                if(snap && is< numprocs && cpuid==0 && it%jt==0)
                    sf_floatwrite(u1[0], padnx*padnz, snapshot);
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(permm[0], 0, nz*nx*sizeof(float));
            
            v1=is*ds_v;
            v2=v1+nx;
            for(ix=v1; ix<v2; ix++){
                for(iz=0; iz<nz; iz++){
                    permm[ix-v1][iz]+=mm[ix][iz]/(sou2[ix-v1][iz]+FLT_EPSILON);
                }
            }
            
            for(it=0; it<nt; it++){
                //sf_warning("Modeling_Receiver: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        u1[ix+padx][iz+padz]+=wave[it][ix][iz]*permm[ix][iz];
                    }
                }
                
                for(ir=0; ir<nr; ir++){
                    rx=ir*dr_v+r0_v;
                    dd[ir][it]+=u1[rx][zr_v];
                }
            } //end of it
            
            fseeko(swap, is*nr*nt*sizeof(float), 0);
            fwrite(dd[0], sizeof(float)*nt, nr, swap);
        }// end of shot
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(cpuid==0){
            for(is=0; is<ns; is++){
                fseeko(swap, is*nr*nt*sizeof(float), 0);
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
    double dt2, idx2, idz2;
 
    float **mm;
    sf_file vel, wavelet;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    sf_init(argc, argv);
    omp_init();
    
    if(cpuid==0) sf_warning("numprocs=%d", numprocs);
    
    if(!sf_getbool("adj", &adj)) adj=true;
    if(!sf_getbool("verb", &verb)) verb=false;
    if(!sf_getbool("snap", &snap)) snap=false;
    
    in=sf_input("input");
    out=sf_output("output");
    vel=sf_input("velocity");
    wavelet=sf_input("wavelet");
    
    if(!sf_histint(vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "o1", &z0)) sf_error("No o1= in velocity");
    
    if(!sf_histint(vel, "n2", &nvx)) sf_error("No n2= in velocity");
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
        
        if(cpuid==0){
        sf_putint(out, "n1", nz);
        sf_putfloat(out, "d1", dz);
        sf_putfloat(out, "o1", z0);
        sf_putint(out, "n2", nvx);
        sf_putfloat(out, "d2", dx);
        sf_putfloat(out, "o2", x0);
        sf_putint(out, "n3", 1);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "m");
        sf_putstring(out, "label2", "Distance");
        sf_putstring(out, "unit2", "m");}
    }else{/* modeling */
        if(!sf_getint("nr", &nr)) sf_error("No nr=");
        if(!sf_getfloat("dr", &dr)) sf_error("No dr=");
        if(!sf_getfloat("r0", &r0)) sf_error("No r0=");
        
        if(!sf_getint("ns", &ns)) sf_error("No ns=");
        if(!sf_getfloat("ds", &ds)) sf_error("No ds=");
        if(!sf_getfloat("s0", &s0)) sf_error("No s0=");
        
        if(cpuid==0){
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
        sf_putstring(out, "label2", "Distance");
        sf_putstring(out, "unit2", "m");
        sf_putstring(out, "label3", "Shot");
        sf_putstring(out, "unit3", "m");}
    }
    
    if(!sf_getfloat("zr", &zr)) zr=0.0;
    if(!sf_getfloat("zs", &zs)) zs=0.0;
    if(!sf_getint("jt", &jt)) jt=100;
    if(!sf_getint("nx", &nx)) sf_error("Need nx");
    if(!sf_getint("padz", &padz)) sf_error("Need padz");
    if(!sf_getint("padx", &padx)) sf_error("Need padx");
    if((s0-x0)/dx*2+1 != nx) sf_error("Need to reset nx");
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    padx0=x0-dx*padx;
    padz0=z0-dz*padz;
    
    dr_v=(dr/dx)+0.5;
    r0_v=(r0-padx0)/dx+0.5;
    zr_v=(zr-padz0)/dz+0.5;
    
    ds_v=(ds/dx)+0.5;
    sx=(s0-padx0)/dx+0.5;
    zs_v=(zs-padz0)/dz+0.5;
    
    mm=sf_floatalloc2(nz, nvx);
    vv=sf_floatalloc2(nz, nvx);
    padvv=sf_floatalloc2(padnz, padnx);
    ww=sf_floatalloc(nt);
    
    sf_floatread(ww, nt, wavelet);
    sf_floatread(vv[0], nz*nvx, vel);
    
    dt2=dt*dt;
    for(ix=0; ix<nvx; ix++)
        for(iz=0; iz<nz; iz++)
            vv[ix][iz]=vv[ix][iz]*vv[ix][iz]*dt2;
    
    if(snap){
        snapshot=sf_output("snapshot");
        
        if(cpuid==0){
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
        sf_putstring(snapshot, "unit1", "m");
        sf_putstring(snapshot, "label2", "Distance");
        sf_putstring(snapshot, "unit2", "m");
        sf_putstring(snapshot, "label3", "Time");
        sf_putstring(snapshot, "unit3", "s");}
    }
    
    if(cpuid==0){
    sf_warning("nx=%d nvx=%d nz=%d nt=%d nr=%d ns=%d", nx, nvx, nz, nt, nr, ns);
    sf_warning("padnx=%d padnz=%d padx0=%.3f padz0=%.3f", padnx, padnz, padx0, padz0);
    sf_warning("dx=%.3f dz=%.3f dt=%.3f dr=%.3f ds=%.3f", dx, dz, dt, dr, ds);
    sf_warning("x0=%.3f z0=%.3f t0=%.3f r0=%.3f s0=%.3f", x0, z0, t0, r0, s0);
    sf_warning("dr_v=%d r0_v=%d zr_v=%d ds_v=%d sx=%d zs_v=%d", dr_v, r0_v, zr_v, ds_v, sx, zs_v);}
    
    idz2=1./(dz*dz);
    idx2=1./(dx*dx);
//    c0=-5.0/2.;
//    c1=4.0/3.;
//    c2=-1.0/12.;
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2*(c11+c12+c21+c22);
    
    if(adj){/* migration */
        prertm2_oper(adj, mm);
        if(cpuid==0) sf_floatwrite(mm[0], nz*nvx, out);
    }else{
        sf_floatread(mm[0], nz*nvx, in);
        prertm2_oper(adj, mm);
    }
    
    sf_fileclose(in); sf_fileclose(out);    
    free(*vv); free(*padvv);
    free(vv); free(padvv); free(ww);
    MPI_Finalize();
    
    exit(0);
}
