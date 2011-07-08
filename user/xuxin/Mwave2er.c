/* 2D acoustic exploding reflector modeling */
/* Equation: P_{tt} = K * div(1/rho* grad(P))
   Initial:  P(x,t=0)=I(x) and P_t(t=0)=0
   Boundary: one-way ABC + sponge
   Receiver: (xr1,zr1) ...
   FD:       4th space and 2nd time 
*/

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
#define DXDX(a,ix,iz,idx2) (C0* a[iz][ix]+\
                            CA*(a[iz][ix-1]+a[iz][ix+1])+\
                            CB*(a[iz][ix-2]+a[iz][ix+2]))*idx2
#define DZDZ(a,ix,iz,idz2) (C0* a[iz][ix]+\
                            CA*(a[iz-1][ix]+a[iz+1][ix])+\
                            CB*(a[iz-2][ix]+a[iz+2][ix]))*idz2
#define C1  0.666666666666666666667 /*  2/3 */
#define C2 -0.083333333333333333333 /* -1/12 */
#define DX(a,ix,iz,idx) (C1*(a[iz][ix+1]-a[iz][ix-1])+\
                         C2*(a[iz][ix+2]-a[iz][ix-2]))*idx
#define DZ(a,ix,iz,idz) (C1*(a[iz+1][ix]-a[iz-1][ix])+\
                         C2*(a[iz+2][ix]-a[iz-2][ix]))*idz

int main(int argc, char* argv[])
{
    sf_file Fwav;    float  *wav=NULL;
    sf_file Fvel;    float **vel=NULL;
    sf_file Fden;    float **den=NULL,**rox=NULL,**roz=NULL;
    sf_file Fimg;    float **img=NULL;
    sf_file Frec;    float  *dd =NULL;
    sf_file Fdat;    float **dat=NULL;
    sf_file Fout;    float **um =NULL,**uo=NULL,**up=NULL,**ua=NULL;
    int nx,ix,nbx,npx,x; float ox,dx,idx,idx2;
    int nz,iz,nbz,npz,z; float oz,dz,idz,idz2;
    int nt,it;           float ot,dt,dt2;
    int nr;              float *xr,*zr;
    bool fsrf;
    int jsnap,nsnap,ompnth,ompath,ompchunk,iz_srf;
    float **vt2=NULL,**tt=NULL;
    wave2dp wave=NULL;
    abcone2dp abc=NULL;
    sponge2dp spo=NULL;
    lint2dp lintr=NULL;

    /* read parameters */
    sf_init(argc,argv);
    Fwav=sf_input("in");    /* [1][nt]         source wavelet */ 
    Fvel=sf_input("vel");   /* [nz][nx]        velocity */
    Fden=sf_input("den");   /* [nz][nx]        density */
    Fimg=sf_input("img");   /* [nz][nx]        reflectivity */
    Frec=sf_input("rec");   /* [2][nr]         receiver coord */
    Fout=sf_output("out");  /* [nsnap][nz][nx] modeled wavefield*/
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
    dim2d_equal(Fimg,"img",nx,ox,dx,nz,oz,dz);
    if (!sf_histint(Frec,"n1",&nr))   sf_error("No n1= in rec");
    if (!sf_getbool("free",&fsrf))    fsrf=false;
    /* free surface */
    if (!sf_getint("jsnap",&jsnap))   jsnap=1;
    /* jump in wavefield */
    if (!sf_getint("nb1",&nbx))       nbx=20;
    /* num of boundary grid in axis 1 */
    if (!sf_getint("nb2",&nbz))       nbz=20;
    /* num of boundary grid in axis 2 */
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
    sf_putstring(Fout,"label2","z");
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
    npx   = wave->npx;
    npz   = wave->npz;
    iz_srf= wave->iz_srf;

    wav=sf_floatalloc(nt);
    vel=sf_floatalloc2(npx,npz);
    vt2=sf_floatalloc2(npx,npz);
    img=sf_floatalloc2(npx,npz);
    den=sf_floatalloc2(npx,npz);
    rox=sf_floatalloc2(npx,npz);
    roz=sf_floatalloc2(npx,npz);
    um =sf_floatalloc2(npx,npz);
    uo =sf_floatalloc2(npx,npz);
    up =sf_floatalloc2(npx,npz);
    ua =sf_floatalloc2(npx,npz);
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
    /* This density term is nonconservative. Staggered grid is a conservative implementation but I got dispered reflection. So far the nonconservative method is ok. */
    for     (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    x=ix+nbx; z=iz+iz_srf;
	    rox[z][x] = DX(den,x,z,idx)/den[z][x];
	    if (fsrf && iz<2) roz[z][x]=0.0f;
	    else roz[z][x] = DZ(den,x,z,idz)/den[z][x];
	}
    }

    /* read reflectivity */
    sf_floatread(tt[0],nx*nz,Fimg);
    expand(tt,img,wave);

    /* initial value */ 
    for     (iz=0; iz<npz; iz++) 
	for (ix=0; ix<npx; ix++) 
	    vt2[iz][ix] = vel[iz][ix]*vel[iz][ix] * dt2;
    
    for (ix=0; ix<npx*npz; ix++) {
	um[0][ix] = 0.0f;
	uo[0][ix] = 0.0f;
	up[0][ix] = 0.0f;
	ua[0][ix] = 0.0f;
    }

    abc=abcone2d_init(vel,wave);
    spo=sponge2d_init(0.00001,wave);

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
	/* laplacian */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(uo,ua,nx,nz,nbx,iz_srf,idx2,idz2,fsrf,rox,roz)
#endif	
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x=ix+nbx;
		z=iz+iz_srf;
		if (fsrf && iz<2)
		    ua[z][x]=0.0f;
		else
		    ua[z][x] = DXDX(uo,x,z,idx2)
                             + DZDZ(uo,x,z,idz2)
                             - rox[z][x]*DX(uo,x,z,idx)
                             - roz[z][x]*DZ(uo,x,z,idz);
	    }
	}
	/* time stepping */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,um,up,vt2,nx,nz,nbx,iz_srf,fsrf)
#endif
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x=ix+nbx; 
		z=iz+iz_srf;
		if (fsrf && iz<2)
		    up[z][x]=0.0f;
		else
		    up[z][x]=2.0f*uo[z][x]-um[z][x]+vt2[z][x]*ua[z][x];
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

/* Error: **glibc detected** 
    free(*um);   free(um);
    free(*uo);   free(uo);
    free(*up);   free(up);
*/
    free(*rox);  free(rox);
    free(*roz);  free(roz);
    free(*ua);   free(ua);
    free(*vt2);  free(vt2);
    free(*tt);   free(tt);
    free(dd);
    free(*dat);  free(dat);
    free(*vel);  free(vel);
    free(*den);  free(den);
    free(*img);  free(img);

    exit(0);
}
