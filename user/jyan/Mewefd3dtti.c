/* 3D elastic time-domain FD modeling */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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

#include "fdutil3d.h"

/* 
 * axes: z(1), x(2), y(3)
 */


#define NOP 4 /* derivative operator half-size */

#define C1 +0.800000   /* +4/5    */
#define C2 -0.200000   /* -1/5    */
#define C3 +0.038095   /* +4/105  */
#define C4 -0.003571   /* -5/280  */
#define Dx(a,ix,iy,iz,s) (C4*(a[iy][ix+4][iz] - a[iy][ix-4][iz]) + \
			  C3*(a[iy][ix+3][iz] - a[iy][ix-3][iz]) +   \
			  C2*(a[iy][ix+2][iz] - a[iy][ix-2][iz]) +	\
			  C1*(a[iy][ix+1][iz] - a[iy][ix-1][iz])  )*s
#define Dy(a,ix,iy,iz,s) (C4*(a[iy+4][ix][iz] - a[iy-4][ix][iz]) + \
			  C3*(a[iy+3][ix][iz] - a[iy-3][ix][iz]) +   \
			  C2*(a[iy+2][ix][iz] - a[iy-2][ix][iz]) +	\
			  C1*(a[iy+1][ix][iz] - a[iy-1][ix][iz])  )*s
#define Dz(a,ix,iy,iz,s) (C4*(a[iy][ix][iz+4] - a[iy][ix][iz-4]) + \
			  C3*(a[iy][ix][iz+3] - a[iy][ix][iz-3]) +   \
			  C2*(a[iy][ix][iz+2] - a[iy][ix][iz-2]) +	\
			  C1*(a[iy][ix][iz+1] - a[iy][ix][iz-1])  )*s
/*------------------------------------------------------------*/


int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,ssou,dabc;
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,ax,ay,az;
    sf_axis as,ar,ac;

    int     nt,nz,nx,ny,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy,idz,idx,idy;
  

    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abcp=NULL,abcs=NULL;
    sponge1d  spo=NULL;

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    float **dd=NULL;           /* data      */

    /*------------------------------------------------------------*/
    float ***tt=NULL;
    float ***ro=NULL;           /* density */

    /* orthorombic stiffness - 9 coefficients */
    /* c11 c12 c13 
       .   c22 c23 
       .   .   c33 
                  c44
                     c55
                        c66 */
    float ***c11=NULL,***c12=NULL,***c13=NULL,***c14=NULL,***c15=NULL,***c16=NULL;
    float ***c22=NULL,***c23=NULL,***c24=NULL,***c25=NULL,***c26=NULL;
    float ***c33=NULL,***c34=NULL,***c35=NULL,***c36=NULL;
    float ***c44=NULL,***c45=NULL,***c46=NULL;
    float ***c55=NULL,***c56=NULL;
    float ***c66=NULL;
  
    float ***vp,***vs;

    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float ***umz,***uoz,***upz,***uaz,***utz; 
    float ***umx,***uox,***upx,***uax,***utx;
    float ***umy,***uoy,***upy,***uay,***uty;

    /* stress/strain tensor */ 
    float ***tzz,***txx,***tyy,***txy,***tyz,***tzx;       
    float    szz,   sxx,   syy,   sxy,   syz,   szx;

    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices */
    lint3d cs,cr;

    /* */
    int nbell;
    
    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***uc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* execution flags */
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"time");     sf_setunit(at,"s"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"space z");  sf_setunit(az,"km");if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"space x");  sf_setunit(ax,"km");if(verb) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"space y");  sf_setunit(ay,"km");if(verb) sf_raxa(ay); /* space y */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"sources");  sf_setunit(as,"km");if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"receivers");sf_setunit(ar,"km");if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    ns = sf_n(as);
    nr = sf_n(ar);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }


    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    fdbell3d_init(nbell);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); sf_setlabel(az,"expanded z");sf_setunit(az,"km");if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); sf_setlabel(ax,"expanded x");sf_setunit(ax,"km");if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); sf_setlabel(ay,"expanded y");sf_setunit(ay,"km");if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

    /* 3D vector components */
    nc=3;
    ac=sf_maxa(nc,0,1);

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);
	if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
	if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);

	dqz=sf_d(az);
	dqx=sf_d(ax);
	dqy=sf_d(ay);

	acz = sf_maxa(nqz,oqz,dqz);  sf_setunit(acz,"km");sf_setlabel(acz,"snapshot z");sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx);  sf_setunit(acx,"km");sf_setlabel(acx,"snapshot x");sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy);  sf_setunit(acy,"km");sf_setlabel(acy,"snapshot y");sf_raxa(acy);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	sf_setlabel(at,"snapshot frames");
	if(verb)  sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,acy,3);
	sf_oaxa(Fwfl,ac, 4);
	sf_oaxa(Fwfl,at, 5);
    }

    /*------------------------------------------------------------*/
    /* source array */
    ww=sf_floatalloc3(ns,nc,nt); 
    sf_floatread(ww[0][0],nt*nc*ns,Fwav);

    /* data array */
    dd=sf_floatalloc2(nr,nc);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */
    
    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 
    
    ro =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    c11=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c12=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c13=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c14=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c15=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c16=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c22=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c23=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c24=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);    
    c25=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c26=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c33=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c34=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c35=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c36=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c44=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c45=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c46=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);    
    c55=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c56=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c66=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);   

 

    /* input density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden);     expand3d(tt,ro ,fdm);

    /* input stiffness */
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c14,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c15,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c16,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c24,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c25,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c26,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c34,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c35,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c36,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c45,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c46,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c56,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,fdm);

    free(**tt); free(*tt); free(tt);

    /*------------------------------------------------------------*/
    if(dabc) {
	/* one-way abc setup   */
	vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	vs = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    vp[iy][ix][iz] = sqrt( c11[iy][ix][iz]/ro[iy][ix][iz] );
		    vs[iy][ix][iz] = sqrt( c13[iy][ix][iz]/ro[iy][ix][iz] );
		}
	    }
	}
	abcp = abcone3d_make(NOP,dt,vp,fsrf,fdm);
	abcs = abcone3d_make(NOP,dt,vs,fsrf,fdm);
	free(**vp); free(*vp); free(vp);
	free(**vs); free(*vs); free(vs);

	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		ro[iy][ix][iz] = dt*dt/ro[iy][ix][iz];
	    }
	}
     }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    umz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uaz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uax=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uay=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    tzz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tzx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		umz[iy][ix][iz]=0; umx[iy][ix][iz]=0; umy[iy][ix][iz]=0;
		uoz[iy][ix][iz]=0; uox[iy][ix][iz]=0; uoy[iy][ix][iz]=0;
		upz[iy][ix][iz]=0; upx[iy][ix][iz]=0; upy[iy][ix][iz]=0;
		uaz[iy][ix][iz]=0; uax[iy][ix][iz]=0; uay[iy][ix][iz]=0;
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");

    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

	/*------------------------------------------------------------*/
	/* from displacement to strain                                */
	/*------------------------------------------------------------*/	
	/* 
	 * exx = Dx(ux)
	 * eyy = Dy(uy)
	 * ezz = Dz(uz)
	 * exy = Dy(ux) + Dx(uy)
	 * eyz = Dz(uy) + Dy(uz)
	 * ezx = Dx(uz) + Dz(ux)
	 */
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,idx,idy,idz)
#endif
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		    
		    txx[iy][ix][iz] = Dx(uox,ix,iy,iz,idx);
		    tyy[iy][ix][iz] = Dy(uoy,ix,iy,iz,idy);
		    tzz[iy][ix][iz] = Dz(uoz,ix,iy,iz,idz);
		    
		    txy[iy][ix][iz] = Dy(uox,ix,iy,iz,idy) + Dx(uoy,ix,iy,iz,idx);
		    tyz[iy][ix][iz] = Dz(uoy,ix,iy,iz,idz) + Dy(uoz,ix,iy,iz,idy);
		    tzx[iy][ix][iz] = Dx(uoz,ix,iy,iz,idx) + Dz(uox,ix,iy,iz,idz);
		}
	    }
	}
	
	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    
		    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
			+ c12[iy][ix][iz] * tyy[iy][ix][iz]
			+ c13[iy][ix][iz] * tzz[iy][ix][iz]
			+ c14[iy][ix][iz] * tyz[iy][ix][iz]
			+ c15[iy][ix][iz] * tzx[iy][ix][iz]
			+ c16[iy][ix][iz] * txy[iy][ix][iz];

		    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
			+ c22[iy][ix][iz] * tyy[iy][ix][iz]
			+ c23[iy][ix][iz] * tzz[iy][ix][iz]
			+ c24[iy][ix][iz] * tyz[iy][ix][iz]
			+ c25[iy][ix][iz] * tzx[iy][ix][iz]
			+ c26[iy][ix][iz] * txy[iy][ix][iz];

		    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
			+ c23[iy][ix][iz] * tyy[iy][ix][iz]
			+ c33[iy][ix][iz] * tzz[iy][ix][iz]
			+ c34[iy][ix][iz] * tyz[iy][ix][iz]
			+ c35[iy][ix][iz] * tzx[iy][ix][iz]
			+ c36[iy][ix][iz] * txy[iy][ix][iz];

		    syz = c14[iy][ix][iz] * txx[iy][ix][iz]
			+ c24[iy][ix][iz] * tyy[iy][ix][iz]
			+ c34[iy][ix][iz] * tzz[iy][ix][iz]
			+ c44[iy][ix][iz] * tyz[iy][ix][iz]
			+ c45[iy][ix][iz] * tzx[iy][ix][iz]
			+ c46[iy][ix][iz] * txy[iy][ix][iz];
		    
		    szx = c15[iy][ix][iz] * txx[iy][ix][iz]
			+ c25[iy][ix][iz] * tyy[iy][ix][iz]
			+ c35[iy][ix][iz] * tzz[iy][ix][iz]
			+ c45[iy][ix][iz] * tyz[iy][ix][iz]
			+ c55[iy][ix][iz] * tzx[iy][ix][iz]
			+ c56[iy][ix][iz] * txy[iy][ix][iz];

		    sxy = c16[iy][ix][iz] * txx[iy][ix][iz]
			+ c26[iy][ix][iz] * tyy[iy][ix][iz]
			+ c36[iy][ix][iz] * tzz[iy][ix][iz]
			+ c46[iy][ix][iz] * tyz[iy][ix][iz]
			+ c56[iy][ix][iz] * tzx[iy][ix][iz]
			+ c66[iy][ix][iz] * txy[iy][ix][iz];

 
		  


		    txx[iy][ix][iz] = sxx;
		    tyy[iy][ix][iz] = syy;
		    tzz[iy][ix][iz] = szz;
		    txy[iy][ix][iz] = sxy;
		    tyz[iy][ix][iz] = syz;
		    tzx[iy][ix][iz] = szx;
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	if(fsrf) {
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx)
#endif
	    for        (iy=0; iy<fdm->nypad; iy++) {
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nb;    iz++) {
			txx[iy][ix][iz]=0;
			tyy[iy][ix][iz]=0;
			tzz[iy][ix][iz]=0;
			txy[iy][ix][iz]=0;
			tyz[iy][ix][iz]=0;
			tzx[iy][ix][iz]=0;
		    }
		}
	    }
	}


	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(ssou) {
	    lint3d_bell(txx,ww[it][0],cs);
	    lint3d_bell(tyy,ww[it][1],cs);
	    lint3d_bell(tzz,ww[it][2],cs);
	}

	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* 
	 * ax = Dx(txx) + Dy(txy) + Dz(txz)
	 * ay = Dx(txy) + Dy(tyy) + Dz(tyz)
	 * az = Dx(txz) + Dy(tyz) + Dz(tzz)
	 */	
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,idx,idy,idz)
#endif
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		    
		    uax[iy][ix][iz] = Dx( txx,ix,iy,iz,idx) + Dy( txy,ix,iy,iz,idy) + Dz( tzx,ix,iy,iz,idz) ;
		    uay[iy][ix][iz] = Dx( txy,ix,iy,iz,idx) + Dy( tyy,ix,iy,iz,idy) + Dz( tyz,ix,iy,iz,idz) ;
		    uaz[iy][ix][iz] = Dx( tzx,ix,iy,iz,idx) + Dy( tyz,ix,iy,iz,idy) + Dz( tzz,ix,iy,iz,idz) ;		    
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(!ssou) {
	    lint3d_bell(uaz,ww[it][0],cs);
	    lint3d_bell(uax,ww[it][1],cs);
	    lint3d_bell(uay,ww[it][2],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,uox,uoy,uoz,umx,umy,umz,upx,upy,upz,uax,uay,uaz,ro)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    upx[iy][ix][iz] = 2*uox[iy][ix][iz] 
			-               umx[iy][ix][iz] 
			+               uax[iy][ix][iz] * ro[iy][ix][iz]; 

		    upy[iy][ix][iz] = 2*uoy[iy][ix][iz] 
			-               umy[iy][ix][iz] 
			+               uay[iy][ix][iz] * ro[iy][ix][iz]; 

		    upz[iy][ix][iz] = 2*uoz[iy][ix][iz] 
			-               umz[iy][ix][iz] 
			+               uaz[iy][ix][iz] * ro[iy][ix][iz]; 
		    
		}
	    }
	}
	/* circulate wavefield arrays */
	utz=umz; uty=umy; utx=umx;
	umz=uoz; umy=uoy; umx=uox;
	uoz=upz; uoy=upy; uox=upx;
	upz=utz; upy=uty; upx=utx;
	
	if(dabc) {
	    /* one-way ABC */
	    abcone3d_apply(uoz,umz,NOP,abcp,fdm);
	    abcone3d_apply(uox,umx,NOP,abcp,fdm);
	    abcone3d_apply(uoy,umy,NOP,abcp,fdm);
	    
	    abcone3d_apply(uoz,umz,NOP,abcs,fdm);
	    abcone3d_apply(uox,umx,NOP,abcs,fdm);
	    abcone3d_apply(uoy,umy,NOP,abcs,fdm);

	    /* sponge ABC */
	    sponge3d_apply(umz,spo,fdm);
	    sponge3d_apply(uoz,spo,fdm);
	    sponge3d_apply(upz,spo,fdm);
	    
	    sponge3d_apply(umx,spo,fdm);
	    sponge3d_apply(uox,spo,fdm);
	    sponge3d_apply(upx,spo,fdm);

	    sponge3d_apply(umy,spo,fdm);
	    sponge3d_apply(uoy,spo,fdm);
	    sponge3d_apply(upy,spo,fdm);
	}	    

	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
	lint3d_extract(uoz,dd[0],cr);
	lint3d_extract(uox,dd[1],cr);
	lint3d_extract(uoy,dd[2],cr);

	if(snap && it%jsnap==0) {
	    cut3d(uoz,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

	    cut3d(uox,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

	    cut3d(uoy,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
	}
	if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);

    }
    if(verb) fprintf(stderr,"\n");    
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);

    free(**ro);  free(*ro);  free(ro);
    free(**c11); free(*c11); free(c11);
    free(**c12); free(*c12); free(c12);
    free(**c13); free(*c13); free(c13);
    free(**c14); free(*c14); free(c14);
    free(**c15); free(*c15); free(c15);
    free(**c16); free(*c16); free(c16);
    free(**c22); free(*c22); free(c22);
    free(**c23); free(*c23); free(c23);
    free(**c24); free(*c24); free(c24);
    free(**c25); free(*c25); free(c25);
    free(**c26); free(*c26); free(c26);
    free(**c33); free(*c33); free(c33);
    free(**c34); free(*c34); free(c34);
    free(**c35); free(*c35); free(c35);
    free(**c36); free(*c36); free(c36);
    free(**c44); free(*c44); free(c44);
    free(**c45); free(*c45); free(c45);
    free(**c46); free(*c46); free(c46);
    free(**c55); free(*c55); free(c55);
    free(**c56); free(*c56); free(c56);
    free(**c66); free(*c66); free(c66);




    free(**umz); free(*umz); free(umz);
    free(**uoz); free(*uoz); free(uoz);
    free(**upz); free(*upz); free(upz);
    free(**uaz); free(*uaz); free(uaz);

    free(**umx); free(*umx); free(umx);
    free(**uox); free(*uox); free(uox);
    free(**upx); free(*upx); free(upx);
    free(**uax); free(*uax); free(uax);

    free(**umy); free(*umy); free(umy);
    free(**uoy); free(*uoy); free(uoy);
    free(**upy); free(*upy); free(upy);
    free(**uay); free(*uay); free(uay);

    free(**tzz); free(*tzz); free(tzz);
    free(**txx); free(*txx); free(txx);
    free(**tyy); free(*tyy); free(tyy);
    free(**txy); free(*txy); free(txy);
    free(**tyz); free(*tyz); free(tyz);
    free(**tzx); free(*tzx); free(tzx);

    free(**uc);  free(*uc);  free(uc);    

    exit (0);
}

