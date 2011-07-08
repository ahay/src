/* 2D pseudo-depth exploding reflector modeling */
/* Equation: P_{tt} = v^2 *(P_{xx} + 2*s*P_{xz} + s^2*P_{zz}) + v0^2*P_{zz}
   Initial:  P(t=0)=I(x) and P_t(t=0)=0
   Boundary: one-way ABC + sponge
   Receiver: (xr1,zr1) ...
   FD:       4th space and 2nd time 
*/
/* density is currently ignored */

/* Copyright (C) 2011 KAUST */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "wave.h"

/* FD stencils */
#define C0 -2.500000000000000000000 /* -30/12 */
#define CA  1.333333333333333333333 /*  16/12 */
#define CB -0.083333333333333333333 /*  -1/12 */
#define DXDX(a,ix,iz,idx2) (C0*a[iz][ix]+\
                            CA*(a[iz][ix-1]+a[iz][ix+1])+\
                            CB*(a[iz][ix-2]+a[iz][ix+2]))*idx2
#define DZDZ(a,ix,iz,idz2) (C0*a[iz][ix]+\
                            CA*(a[iz-1][ix]+a[iz+1][ix])+\
                            CB*(a[iz-2][ix]+a[iz+2][ix]))*idz2
#define CX  0.250000000000000000000
#define DXDZ(a,ix,iz,idz,idx) (a[iz-1][ix-1]-a[iz-1][ix+1]-a[iz+1][ix-1]+a[iz+1][ix+1])*idx*idz*CX
#define C1  0.666666666666666666667 /*  2/3 */
#define C2 -0.083333333333333333333 /* -1/12 */
#define DX(a,ix,iz,idx) (C1*(a[iz][ix+1]-a[iz][ix-1])+\
                         C2*(a[iz][ix+2]-a[iz][ix-2]))*idx
#define DZ(a,ix,iz,idz) (C1*(a[iz+1][ix]-a[iz-1][ix])+\
                         C2*(a[iz+2][ix]-a[iz-2][ix]))*idz

int main(int argc, char* argv[])
{
    sf_file Fwav;         float  *wav=NULL;
    sf_file Fvel;         float **vel=NULL;
    sf_file Fden;         float **den=NULL,**rox=NULL,**roz=NULL;
    sf_file Fsig;         float **sig=NULL;
    sf_file Fimg;         float **img=NULL;
    sf_file Frec;         float  *dd =NULL;
    sf_file Fdat;         float **dat=NULL;
    sf_file Fout;         float **um =NULL,**uo=NULL,**up=NULL,**ua=NULL;
    int nx,ix,nbx,npx,x;  float ox,dx,idx,idx2;
    int nz,iz,nbz,npz,z;  float oz,dz,idz,idz2;
    int nt,it;            float ot,dt,dt2;
    int nr;               float *xr,*zr;
    bool fsrf;
    int jsnap,nsnap,ompnth,ompath,ompchunk,iz_srf;
    float **tt=NULL,**cxx=NULL,**cxz=NULL,**czz=NULL,v0,v2;
    wave2dp wave=NULL;
    abcone2dp abc=NULL;
    sponge2dp spo=NULL;
    lint2dp lintr=NULL;

    /* read parameters */
    sf_init(argc,argv);
    Fwav=sf_input("in");    /* [1][nt]        source wavelet (inv=n) */ 
    Fvel=sf_input("vel");   /* [nz][nx]        velocity */
    Fden=sf_input("den");   /* [nz][nx]        density */
    Fsig=sf_input("sig");   /* [nz][nx]        sigma */
    Fimg=sf_input("img");   /* [nz][nx]        reflectivity */
    Frec=sf_input("rec");   /* [2][nr]         receiver coord */
    Fout=sf_output("out");  /* [nsnap][nz][nx] modeled wavefield */
    Fdat=sf_output("dat");  /* [nt][nr]        modeled data */
    if (SF_FLOAT != sf_gettype(Fwav)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fvel)) sf_error("Need float vel");
    if (SF_FLOAT != sf_gettype(Fden)) sf_error("Need float den");
    if (!sf_histint(Fwav,"n1",&nt))   sf_error("No n1= in input");
    if (!sf_histfloat(Fwav,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(Fwav,"d1",&dt)) sf_error("No d1= in input"); dt2=dt*dt;
    if (!sf_histint(Fvel,"n1",&nx))   sf_error("No n1= in vel");
    if (!sf_histfloat(Fvel,"o1",&ox)) sf_error("No o1= in vel");
    if (!sf_histfloat(Fvel,"d1",&dx)) sf_error("No d1= in vel"); idx=1.0f/dx; idx2=idx*idx;
    if (!sf_histint(Fvel,"n2",&nz))   sf_error("No n2= in vel");
    if (!sf_histfloat(Fvel,"o2",&oz)) sf_error("No o2= in vel");
    if (!sf_histfloat(Fvel,"d2",&dz)) sf_error("No d2= in vel"); idz=1.0f/dz; idz2=idz*idz;
    dim2d_equal(Fden,"den",nx,ox,dx,nz,oz,dz);
    dim2d_equal(Fden,"den",nx,ox,dx,nz,oz,dz);
    if (!sf_histint(Frec,"n1",&nr))   sf_error("No n1= in rec");
    if (!sf_getbool("free",&fsrf))    fsrf=false;
    /* free surface */
    if (!sf_getint("jsnap",&jsnap))   jsnap=1;
    /* jump in wavefield */
    if (!sf_getint("nb1",&nbx))       nbx=20;
    /* num of boundary grid in axis 1 */
    if (!sf_getint("nb2",&nbz))       nbz=20;
    /* num of boundary grid in axis 2 */
    if (!sf_getfloat("v0",&v0))       v0=1.0f;
    /* stretching velocity */
    if (!sf_getint("ompnth",&ompnth)) ompnth=1;
    /* OpenMP num of threads */
    if (!sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP chunksize */
#ifdef _OPENMP
#pragma omp parallel
    ompath=omp_get_num_threads();
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    /* output header */
    for (nsnap=0,it=0; it<nt; it++)
	if (!(it % jsnap)) nsnap++;
    sf_putint(Fout,"n1",nx);
    sf_putfloat(Fout,"o1",ox);
    sf_putfloat(Fout,"d1",dx);
    sf_putstring(Fout,"label1","x");
    sf_putstring(Fout,"unit1","");
    sf_putint(Fout,"n2",nz);
    sf_putfloat(Fout,"o2",oz);
    sf_putfloat(Fout,"d2",dz);
    sf_putstring(Fout,"label2","z*");
    sf_putstring(Fout,"unit2","");
    sf_putint(Fout,"n3",nsnap);
    sf_putfloat(Fout,"o3",ot);
    sf_putfloat(Fout,"d3",dt*jsnap);
    sf_putstring(Fout,"label3","t");
    sf_putstring(Fout,"unit3","");

    sf_putint(Fdat,"n1",nr);
    sf_putstring(Fdat,"label1","r");
    sf_putstring(Fdat,"unit1","");
    sf_putint(Fdat,"n2",nt);
    sf_putfloat(Fdat,"o2",ot);
    sf_putfloat(Fdat,"d2",dt);
    sf_putstring(Fdat,"label2","t");
    sf_putstring(Fdat,"unit2","");

    wave=wave2d_init(nx,ox,dx,nbx,nz,oz,dz,nbz,nt,dt,0,nr,fsrf);    
    wave2d_info(wave);
    npx=wave->npx;
    npz=wave->npz;
    iz_srf=wave->iz_srf;

    wav=sf_floatalloc(nt);
    vel=sf_floatalloc2(npx,npz);
    den=sf_floatalloc2(npx,npz);
    rox=sf_floatalloc2(npx,npz);
    roz=sf_floatalloc2(npx,npz);
    sig=sf_floatalloc2(npx,npz);
    img=sf_floatalloc2(npx,npz);
    um =sf_floatalloc2(npx,npz);
    uo =sf_floatalloc2(npx,npz);
    up =sf_floatalloc2(npx,npz);
    ua =sf_floatalloc2(npx,npz);
    cxx=sf_floatalloc2(npx,npz); /* v^2 */
    cxz=sf_floatalloc2(npx,npz); /* 2*sig*v^2 */
    czz=sf_floatalloc2(npx,npz); /* v^2*sig^2+v0^2 */
    tt =sf_floatalloc2(nx ,nz );
    dat=sf_floatalloc2(nr ,nt );
    xr =sf_floatalloc(nr);
    zr =sf_floatalloc(nr);
    dd =sf_floatalloc(nr);

    /* read receiver coord */
    sf_floatread(xr,nr,Frec);
    sf_floatread(zr,nr,Frec);
    lintr = lint2d_init(xr,zr,nr,wave);

    /* read wavelet */
    sf_floatread(wav,nt,Fwav);

    /* read velocity */
    sf_floatread(tt[0],nx*nz,Fvel);
    expand(tt,vel,wave);

    /* read density */
    sf_floatread(tt[0],nx*nz,Fden);
    expand(tt,den,wave);
    for     (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    x=ix+nbx; z=iz+iz_srf;
	    rox[z][x] = DX(den,x,z,idx)/den[z][x];
	    if (fsrf && iz<2) roz[z][x]=0.0f;
	    else roz[z][x] = DZ(den,x,z,idz)/den[z][x];
	}
    }

    /* read sigma */
    sf_floatread(tt[0],nx*nz,Fsig);
    expand(tt,sig,wave);

    /* read reflectivity */
    sf_floatread(tt[0],nx*nz,Fimg);
    expand(tt,img,wave);

    /* initial value */ 
    for     (iz=0; iz<npz; iz++) {
	for (ix=0; ix<npx; ix++) {
	    v2= vel[iz][ix]*vel[iz][ix];
	    cxx[iz][ix] = v2;
	    cxz[iz][ix] = v2*sig[iz][ix]*2.0f;
	    czz[iz][ix] = v2*sig[iz][ix]*sig[iz][ix]+v0*v0;
	}
    }
    for (ix=0; ix<npx*npz; ix++) {
	    um[0][ix] = 0.0f;
	    uo[0][ix] = 0.0f;
	    up[0][ix] = 0.0f;
	    ua[0][ix] = 0.0f;
    }

    abc=abcone2d_tau_init(vel,v0,wave);
    spo=sponge2d_init(0.00005,wave);

    /* MAIN */
    for (it=0; it<nt; it++) {
	fprintf(stderr,"\b\b\b\b\b%d",it);
	/* write wavefield */
	if (!(it % jsnap)) {
	    for(iz=iz_srf; iz<npz-nbz; iz++)
		sf_floatwrite(&uo[iz][nbx],nx,Fout);
	}
	/* extract data */
	lint2d_extract(uo,dd,lintr);
	sf_floatwrite(dd,nr,Fdat);
	/* inject source */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(uo,nx,nz,nbx,iz_srf,fsrf,rox,roz,wav,img)
#endif	
	for (iz=0; iz<nz; iz++)
	    for (ix=0; ix<nx; ix++) {
		x=ix+nbx;
		z=iz+iz_srf;
		if (fsrf && iz<2) break;
		else
		    uo[z][x] += wav[it]*img[z][x];
	    }
	/* v^2*laplacian */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,cxx,cxz,czz,nx,nz,nbx,iz_srf,idx2,idz2,fsrf)
#endif
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x=ix+nbx;
		z=iz+iz_srf;
		if (fsrf && iz<2)
		    ua[z][x]=0.0f;
		else
		    ua[z][x] = DXDX(uo,x,z,idx2)*cxx[z][x] 
                	     + DXDZ(uo,x,z,idz,idx)*cxz[z][x]
			     + DZDZ(uo,x,z,idz2)*czz[z][x];
	    }
	}
	/* time stepping */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,um,up,dt2,nx,nz,nbx,iz_srf)
#endif
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x=ix+nbx;
		z=iz+iz_srf;
		if (fsrf && iz<2)
		    up[z][x]=0.0f;
		else 
		    up[z][x]=uo[z][x]*2.0f-um[z][x]+ua[z][x]*dt2;
	    }
	}
	/* absorbing BC */
	abcone2d_absorb(up,uo,wave,abc);
	/* sponge */
	sponge2d_sponge(up,wave,spo);

	/* circulate pointer */
	um=uo;
	uo=up;
	up=um;
    } /* it */
    fprintf(stderr,"\n");

    /* clear */
    sf_fileclose(Fwav);
    sf_fileclose(Fvel);
    sf_fileclose(Fden);
    sf_fileclose(Fimg);
    sf_fileclose(Frec);
    sf_fileclose(Fdat);
    sf_fileclose(Fout);

/*
    free(*um);   free(um);
    free(*uo);   free(uo);
    free(*up);   free(up);
*/
    free(*rox);  free(rox);
    free(*roz);  free(roz);
    free(*ua);   free(ua);
    free(*cxx);  free(cxx);
    free(*cxz);  free(cxz);
    free(*czz);  free(czz);
    free(*tt);   free(tt);
    free(*dat);  free(dat);
    free(*vel);  free(vel);
    free(*den);  free(den);
    free(*sig);  free(sig);
    free(*img);  free(img);
    free(dd);

    exit(0);
}
