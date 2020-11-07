/* 3D elastic recursive integral time extrapolation of decomposed wave modes using KISS-FFT
   sou wavelet  (nx,ny,nc,nt)
   rec data     (nx,ny,nc,nt)
   sou geometry (nc,nx,ny)
   rec geometry (nc,nx,ny)
*/
/*
  Copyright (C) 2016 The University of Texas at Austin
  
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

  The main framework is following Paul Sava's sfewefd3d program.
  This program uses pseudo-spectral method to calculate spatial derivatives.
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ksutil.h"

int main(int argc, char* argv[])
{
    /*------------------------------------------------------------*/
    /* Execution control, I/O files and geometry                  */
    /*------------------------------------------------------------*/
    bool verb,fsrf,snap,back,esou; /* execution flags */
    int  jsnap,ntsnap,jdata; /* jump along axes */
    int  shft; /* time shift for wavefield matching in RTM */

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Frkp=NULL; /* app. rank */
    sf_file Fltp=NULL; /* left mat  */
    sf_file Frtp=NULL; /* right mat */
    sf_file Fwfp=NULL; /* wavefield */
    sf_file Frks=NULL; /* app. rank */
    sf_file Flts=NULL; /* left mat  */
    sf_file Frts=NULL; /* right mat */
    sf_file Fwfs=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,ax,ay,az; /* time, x, y, z */ 
    sf_axis asx,asy,arx,ary,ac;    /* sou, rec-x, rec-y, component */ 

    /* dimension, index and interval */
    int     nt,nz,nx,ny,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy;
    int     nxyz, nk;

    /* FDM and KSP structure */ //!!!JS
    fdm3d    fdm=NULL;
    dft3d    dft=NULL;
    clr3d    clr_p=NULL, clr_s=NULL;

    /* I/O arrays for sou & rec */
    sf_complex***ww=NULL;    /* wavelet   */
    pt3d        *ss=NULL;    /* sources   */
    pt3d        *rr=NULL;    /* receivers */
    sf_complex **dd=NULL;    /* data      */

    /*------------------------------------------------------------*/
    /* displacement: uo = U @ t; up = U @ t+1                     */
    /*------------------------------------------------------------*/
    sf_complex ***uox, ***uoy, ***uoz, **uo;
    sf_complex ***uox_p, ***uoy_p, ***uoz_p, **uo_p;
    sf_complex ***uox_s, ***uoy_s, ***uoz_s, **uo_s;
    /*sf_complex ***upx, ***upy, ***upz, **up;*/

    /*------------------------------------------------------------*/
    /* lowrank decomposition arrays                               */
    /*------------------------------------------------------------*/
    int ntmp, *n2s_p, *n2s_s;
    sf_complex **lt_p, **rt_p, **lt_s, **rt_s;

    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices                       */
    /*------------------------------------------------------------*/
    lint3d cs,cr; /* for injecting source and extracting data */

    /* Gaussian bell */
    int nbell;
    
    /*------------------------------------------------------------*/
    /* wavefield cut params                                       */
    /*------------------------------------------------------------*/
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    sf_complex***uc=NULL; /* tmp array for output wavefield snaps */

    /*------------------------------------------------------------*/
    /* init RSF                                                   */
    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters                                             */
    /*------------------------------------------------------------*/
#ifdef _OPENMP
    omp_init();
#endif

    /*------------------------------------------------------------*/
    /* read execution flags                                       */
    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("back",&back)) back=false; /* backward extrapolation flag (for rtm) */
    if(! sf_getbool("esou",&esou)) esou=false; /* explosive force source */

    /*------------------------------------------------------------*/
    /* I/O files                                                  */
    /*------------------------------------------------------------*/
    Fwav = sf_input ("in" ); /* wavelet   */
    Fdat = sf_output("out"); /* data      */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fccc = sf_input ("ccc"); /* stiffness */
    Frkp = sf_input ("rkp"); /* app. rank */
    Fltp = sf_input ("ltp"); /* left mat  */
    Frtp = sf_input ("rtp"); /* right mat */
    Fwfp = sf_output("wfp"); /* wavefield */
    Frks = sf_input ("rks"); /* app. rank */
    Flts = sf_input ("lts"); /* left mat  */
    Frts = sf_input ("rts"); /* right mat */
    Fwfs = sf_output("wfs"); /* wavefield */

    /*------------------------------------------------------------*/
    /* axes                                                       */
    /*------------------------------------------------------------*/
    at = sf_iaxa(Fwav,4); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space y */

    asx = sf_iaxa(Fsou,2); sf_setlabel(asx,"sx"); if(verb) sf_raxa(asx); /* sources x */
    asy = sf_iaxa(Fsou,3); sf_setlabel(asy,"sy"); if(verb) sf_raxa(asy); /* sources y */
    arx = sf_iaxa(Frec,2); sf_setlabel(arx,"rx"); if(verb) sf_raxa(arx); /* receivers x */
    ary = sf_iaxa(Frec,3); sf_setlabel(ary,"ry"); if(verb) sf_raxa(ary); /* receivers y */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    ns = sf_n(asx)*sf_n(asy);
    nr = sf_n(arx)*sf_n(ary);

    /*------------------------------------------------------------*/
    /* other execution parameters                                 */
    /*------------------------------------------------------------*/
    if(! sf_getint("nbell",&nbell)) nbell=NOP;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }
    if(back) {
        shft = (nt-1)%jsnap;
        sf_warning("For backward extrapolation, make sure nbell(%d)=0",nbell);
    } else shft = 0;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC                     */
    /*------------------------------------------------------------*/
    if( !sf_getint("nb",&nb)) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    if(nbell) fdbell3d_init(nbell);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);

    /*------------------------------------------------------------*/
    /* 3D vector components                                       */
    /*------------------------------------------------------------*/
    nc=3;
    ac=sf_maxa(nc,0,1); /* output 3 cartesian components */

    /*------------------------------------------------------------*/
    /* setup output data header                                   */
    /*------------------------------------------------------------*/
    sf_settype(Fdat,SF_COMPLEX);
    sf_oaxa(Fdat,arx,1);
    sf_oaxa(Fdat,ary,2);
    sf_oaxa(Fdat,ac,3);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,4);

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

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); sf_raxa(acy);

	uc=sf_complexalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0; /* ntsnap = it/jsnap+1; */
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

        sf_settype(Fwfp,SF_COMPLEX);
	sf_oaxa(Fwfp,acz,1);
	sf_oaxa(Fwfp,acx,2);
	sf_oaxa(Fwfp,acy,3);
	sf_oaxa(Fwfp,ac, 4);
	sf_oaxa(Fwfp,at, 5);

        sf_settype(Fwfs,SF_COMPLEX);
	sf_oaxa(Fwfs,acz,1);
	sf_oaxa(Fwfs,acx,2);
	sf_oaxa(Fwfs,acy,3);
	sf_oaxa(Fwfs,ac, 4);
	sf_oaxa(Fwfs,at, 5);
    }

    /*------------------------------------------------------------*/
    /* source and data array                                      */
    /*------------------------------------------------------------*/
    ww=sf_complexalloc3(ns,nc,nt); /* Fast axis: n_sou > n_comp > n_time */
    sf_complexread(ww[0][0],nt*nc*ns,Fwav);

    dd=sf_complexalloc2(nr,nc);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates                          */
    /*------------------------------------------------------------*/
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    /* calculate 3d linear interpolation coef for sou & rec */
    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* allocate and initialize wavefield arrays                   */
    /*------------------------------------------------------------*/
    /* z-component */
    uoz=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz_p=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz_s=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /*upz=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);*/

    /* x-component */
    uox=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox_p=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox_s=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /*upx=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);*/

    /* y-component */
    uoy=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy_p=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy_s=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /*upy=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);*/

    /* wavefield vector */
    uo = (sf_complex**) sf_alloc(3,sizeof(sf_complex*));
    uo[0] = uox[0][0]; uo[1] = uoy[0][0]; uo[2] = uoz[0][0];
    uo_p = (sf_complex**) sf_alloc(3,sizeof(sf_complex*));
    uo_p[0] = uox_p[0][0]; uo_p[1] = uoy_p[0][0]; uo_p[2] = uoz_p[0][0];
    uo_s = (sf_complex**) sf_alloc(3,sizeof(sf_complex*));
    uo_s[0] = uox_s[0][0]; uo_s[1] = uoy_s[0][0]; uo_s[2] = uoz_s[0][0];
    /*up = (sf_complex**) sf_alloc(3,sizeof(sf_complex*));*/
    /*up[0] = upx[0][0]; up[1] = upy[0][0]; up[2] = upz[0][0];*/

    /* initialize fft and lrk */
    dft = dft3d_init(1,false,false,fdm);
    nxyz= fdm->nypad*fdm->nxpad*fdm->nzpad;
    nk  = dft->nky  *dft->nkx  *dft->nkz;

    /*------------------------------------------------------------*/ 
    /* allocation I/O arrays                                      */
    /*------------------------------------------------------------*/ 
    n2s_p = sf_intalloc(9);
    sf_intread(n2s_p,9,Frkp);
    clr_p = clr3d_make2(n2s_p,fdm);

    n2s_s = sf_intalloc(9);
    sf_intread(n2s_s,9,Frks);
    clr_s = clr3d_make2(n2s_s,fdm);

    if (clr_p->n2_max > clr_s->n2_max) clr3d_init(fdm,dft,clr_p);
    else clr3d_init(fdm,dft,clr_s);

    /* check the dimension */
    if (!sf_histint(Fltp,"n1",&ntmp) || ntmp != nxyz)          sf_error("Need n1=%d in left",nxyz);
    if (!sf_histint(Fltp,"n2",&ntmp) || ntmp != clr_p->n2_sum) sf_error("Need n2=%d in left",clr_p->n2_sum);
    if (!sf_histint(Frtp,"n1",&ntmp) || ntmp != nk)            sf_error("Need n1=%d in right",nk);
    if (!sf_histint(Frtp,"n2",&ntmp) || ntmp != clr_p->n2_sum) sf_error("Need n2=%d in right",clr_p->n2_sum);
  
    lt_p = sf_complexalloc2(nxyz,clr_p->n2_sum); 
    rt_p = sf_complexalloc2(nk  ,clr_p->n2_sum); 
    sf_complexread(lt_p[0],nxyz*clr_p->n2_sum,Fltp);
    sf_complexread(rt_p[0],nk  *clr_p->n2_sum,Frtp);

    if (!sf_histint(Flts,"n1",&ntmp) || ntmp != nxyz)          sf_error("Need n1=%d in left",nxyz);
    if (!sf_histint(Flts,"n2",&ntmp) || ntmp != clr_s->n2_sum) sf_error("Need n2=%d in left",clr_s->n2_sum);
    if (!sf_histint(Frts,"n1",&ntmp) || ntmp != nk)            sf_error("Need n1=%d in right",nk);
    if (!sf_histint(Frts,"n2",&ntmp) || ntmp != clr_s->n2_sum) sf_error("Need n2=%d in right",clr_s->n2_sum);

    lt_s = sf_complexalloc2(nxyz,clr_s->n2_sum); 
    rt_s = sf_complexalloc2(nk  ,clr_s->n2_sum); 
    sf_complexread(lt_s[0],nxyz*clr_s->n2_sum,Flts);
    sf_complexread(rt_s[0],nk  *clr_s->n2_sum,Frts);

    /* initialize to zero */
#ifdef _OPENMP
#pragma omp parallel for              \
    schedule(dynamic,1)               \
    private(iy,ix,iz)                 \
    shared(fdm,uoz,uox,uoy,uoz_p,uox_p,uoy_p,uoz_s,uox_s,uoy_s)
#endif
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		uoz  [iy][ix][iz]=sf_cmplx(0.,0.); uox  [iy][ix][iz]=sf_cmplx(0.,0.); uoy  [iy][ix][iz]=sf_cmplx(0.,0.);
		uoz_p[iy][ix][iz]=sf_cmplx(0.,0.); uox_p[iy][ix][iz]=sf_cmplx(0.,0.); uoy_p[iy][ix][iz]=sf_cmplx(0.,0.);
		uoz_s[iy][ix][iz]=sf_cmplx(0.,0.); uox_s[iy][ix][iz]=sf_cmplx(0.,0.); uoy_s[iy][ix][iz]=sf_cmplx(0.,0.);
		/*upz[iy][ix][iz]=sf_cmplx(0.,0.); upx[iy][ix][iz]=sf_cmplx(0.,0.); upy[iy][ix][iz]=sf_cmplx(0.,0.);*/
	    }
	}
    }

    /*------------------------------------------------------------*/ 
    /*------------------------ MAIN LOOP -------------------------*/ 
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
        if(verb) sf_warning("it=%d/%d;",it,nt); /*fprintf(stderr,"\b\b\b\b\b%d",it);*/

	/*------------------------------------------------------------*/
	/* apply lowrank matrix to wavefield vector                   */
        /*------------------------------------------------------------*/
        clr3d_apply(uo_p, uo, lt_p, rt_p, fdm, dft, clr_p);
        clr3d_apply(uo_s, uo, lt_s, rt_s, fdm, dft, clr_s);

	/*------------------------------------------------------------*/
        /* combine P and S wave modes                                 */
        /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for              \
        schedule(dynamic,1)               \
        private(iy,ix,iz)                 \
        shared(fdm,uoz,uox,uoy,uoz_p,uox_p,uoy_p,uoz_s,uox_s,uoy_s)
#endif
        for        (iy=0; iy<fdm->nypad; iy++) {
            for    (ix=0; ix<fdm->nxpad; ix++) {
                for(iz=0; iz<fdm->nzpad; iz++) {
                    uoz[iy][ix][iz] = uoz_p[iy][ix][iz] + uoz_s[iy][ix][iz];
                    uox[iy][ix][iz] = uox_p[iy][ix][iz] + uox_s[iy][ix][iz];
                    uoy[iy][ix][iz] = uoy_p[iy][ix][iz] + uoy_s[iy][ix][iz];
                    /*upz[iy][ix][iz]=sf_cmplx(0.,0.); upx[iy][ix][iz]=sf_cmplx(0.,0.); upy[iy][ix][iz]=sf_cmplx(0.,0.);*/
                }
            }
        }

	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	if(fsrf) { /* need to do something here */ }

        /*------------------------------------------------------------*/
	/* inject displacement source                                 */
	/*------------------------------------------------------------*/
        if(esou) {
            /* exploding force source */
            lint3d_expl_complex(uoz,uox,uoy,ww[it],cs);
        } else {
            if(nbell) {
                lint3d_bell_complex(uoz,ww[it][0],cs);
                lint3d_bell_complex(uox,ww[it][1],cs);
                lint3d_bell_complex(uoy,ww[it][2],cs);
            } else {
                lint3d_inject_complex(uoz,ww[it][0],cs);
                lint3d_inject_complex(uox,ww[it][1],cs);
                lint3d_inject_complex(uoy,ww[it][2],cs);
            }
        }

	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
        if(snap && (it-shft)%jsnap==0) {
            /* P wave */
            cut3d_complex(uoz_p,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfp);

            cut3d_complex(uox_p,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfp);

            cut3d_complex(uoy_p,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfp);

            /* S wave */
            cut3d_complex(uoz_s,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfs);

            cut3d_complex(uox_s,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfs);

            cut3d_complex(uoy_s,uc,fdm,acz,acx,acy);
            sf_complexwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfs);
        }

        lint3d_extract_complex(uoz,dd[0],cr);
        lint3d_extract_complex(uox,dd[1],cr);
        lint3d_extract_complex(uoy,dd[2],cr);
        if(it%jdata==0) sf_complexwrite(dd[0],nr*nc,Fdat);

    }
    if(verb) sf_warning(".");
    if(verb) fprintf(stderr,"\n");    
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    free(dft);
    dft3d_finalize();
    free(clr_p); free(clr_s);
    clr3d_finalize();

    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);

    free(n2s_p); free(n2s_s);
    free(*lt_p); free(lt_p); free(*rt_p); free(rt_p);
    free(*lt_s); free(lt_s); free(*rt_s); free(rt_s);

    free(**uoz  ); free(*uoz  ); free(uoz)  ;
    free(**uoz_p); free(*uoz_p); free(uoz_p);
    free(**uoz_s); free(*uoz_s); free(uoz_s);
    /*free(**upz); free(*upz); free(upz);*/
    free(uo); free(uo_p); free(uo_s);
    /*free(up);*/

    if (snap) {
       free(**uc);  free(*uc);  free(uc);    
    }

    /*------------------------------------------------------------*/


    exit (0);
}

