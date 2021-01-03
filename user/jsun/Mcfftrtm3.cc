// Lowrank decomposition for 3-D isotropic wave propagation. 
//   Copyright (C) 2010 University of Texas at Austin
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <rsf.hh>
#include "lowrank.hh"
/* head files aumatically produced from C programs */
extern "C"{
#include "rtmutil.h"
#include "revolve.h"
}

//using namespace std;

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    /*************************************************************/
    /* controlling flags */
    iRSF par(0);
    bool verb; par.get("verb",verb,false); // verbosity
    bool migr; par.get("migr",migr,false); // adjoint(migration) flag
    bool roll; par.get("roll",roll,false); // rolling v.s. fixed-spread acquisition geometry
    bool dabc; par.get("dabc",dabc,false); // absorbing boundary
    bool lrdc; par.get("lrdc",lrdc,false); // only perform LowRand DeComposition     
    int it,iy,ix,iz,i;

    /*************************************************************/
    /* spatial and time axis */
    sf_file Fvel = sf_input("vel"); // velocity
    sf_file Fwav = sf_input("wav"); // source wavelet

    sf_warning("Input model dimensions...");
    sf_axis az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth z */
    sf_axis ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space x */
    sf_axis ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space y */
    sf_axis at = sf_iaxa(Fwav,1); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time  t */

    int nz = sf_n(az); float dz = sf_d(az); float oz = sf_o(az);
    int nx = sf_n(ax); float dx = sf_d(ax); float ox = sf_o(ax);
    int ny = sf_n(ay); float dy = sf_d(ay); float oy = sf_o(ay);
    int nt = sf_n(at); float dt = sf_d(at); /*float ot = sf_o(at);*/

    /* source */
    sf_complex *ww = sf_complexalloc(nt); /* Fast axis: n_time > n_shot */
    sf_complexread(ww,nt,Fwav);

    /* space domain prep */
    int nb; par.get("nb",nb); // abc width

    fdm3d fdm = fdutil3d_init(verb,false,az,ax,ay,nb,1);
    int nm = fdm->nzpad*fdm->nxpad*fdm->nypad;

    float ***tt  = sf_floatalloc3(nz,nx,ny); 
    float ***vel = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    sf_floatread(tt[0][0],nz*nx*ny,Fvel);
    expand3d(tt,vel,fdm);

    /* set up abc */
    sponge spo=NULL;
    if (dabc) {
	spo = sponge_make(fdm->nb);
    }

    /* wavenumber domain prep */
    dft3d dft = dft3d_init(1,false,false,fdm);
    int nn = dft->nkz*dft->nkx*dft->nky;

    /*************************************************************/
    /* lowrank prep */
    sf_complex **lt, **rt;
    int nr;
    if ( !lrdc && NULL!=sf_getstring("left") && NULL!=sf_getstring("right") ) { /* read in lt and rt */

        sf_file Flft = sf_input("left");
        sf_file Frht = sf_input("right");
        sf_axis am  = sf_iaxa(Flft,1);
        sf_axis arl = sf_iaxa(Flft,2);
        sf_axis arr = sf_iaxa(Frht,1);
        sf_axis an  = sf_iaxa(Frht,2);

        if ( (sf_n(am)!=nm) || (sf_n(an)!=nn) || (sf_n(arl)!=sf_n(arr)) ) sf_error("Rank unmatched!");
        nr = sf_n(arl);

        lt = sf_complexalloc2(nm,nr); sf_complexread(lt[0],nm*nr,Flft);
        rt = sf_complexalloc2(nr,nn); sf_complexread(rt[0],nr*nn,Frht);

    } else { /* perform lowrank decomposition on-the-fly */

        int jump; par.get("jump",jump,1); // subsampling rate for lowrank decomposition
        int seed; par.get("seed",seed,time(NULL)); // seed for random number generator
        float eps; par.get("eps",eps,1.e-4); // tolerance
        int npk; par.get("npk",npk,20); // maximum rank

        lowrank_init(jump, seed, npk, eps, dt, vel[0][0], fdm, dft);
        nr = lowrank_rank();
        lt = sf_complexalloc2(nm,nr);
        rt = sf_complexalloc2(nr,nn);
        lowrank_mat(lt,rt);

    }

    if (lrdc) {
        sf_axis am=NULL, an=NULL, ar=NULL, a1=NULL;
        am = sf_maxa(nm,0,1);  if(verb) sf_raxa(am);
        an = sf_maxa(nn,0,1);    if(verb) sf_raxa(an);
        ar = sf_maxa(nr,0,1); if(verb) sf_raxa(ar);
        a1 = sf_maxa(1,0,1);     if(verb) sf_raxa(a1);

        sf_file Flft = sf_output("left");
        sf_settype(Flft,SF_COMPLEX);
        sf_oaxa(Flft,am,1);
        sf_oaxa(Flft,ar,2);
        sf_oaxa(Flft,a1,3);
        sf_complexwrite(lt[0],nm*nr,Flft);

        sf_file Frht = sf_output("right");
        sf_settype(Frht,SF_COMPLEX);
        sf_oaxa(Frht,ar,1);
        sf_oaxa(Frht,an,2);
        sf_oaxa(Frht,a1,3);
        sf_complexwrite(rt[0],nr*nn,Frht);

        return 0;
    }

    /*************************************************************/
    /* imaging parameters */
    /* input/output -- data and image */
    sf_file Fdat=NULL;
    sf_file Fimg=NULL;
    sf_axis adt = NULL, adx = NULL, ady = NULL;
    int rec_dep, rec_ox, rec_oy, rec_jt, rec_jx, rec_jy, rec_nt, rec_nx, rec_ny;
    par.get("rec_dep",rec_dep,0);
    if (migr) {
        Fdat = sf_input("in");
        sf_warning("Input data dimensions...");
        adt = sf_iaxa(Fdat,1); if(verb) sf_raxa(adt); /* time  t */
        adx = sf_iaxa(Fdat,2); if(verb) sf_raxa(adx); /* space x */
        ady = sf_iaxa(Fdat,3); if(verb) sf_raxa(ady); /* space y */
        if ( sf_o(adx)<ox || sf_o(adx)>ox+nx*dx || sf_o(ady)<oy || sf_o(ady)>oy+ny*dy ) sf_error("Data dimension unmatched!");
        rec_nt = sf_n(adt); rec_jt = (int)(sf_d(adt)/dt); /*rec_ot = (int)((sf_o(adt)-ot)/dt);*/
        rec_nx = sf_n(adx); rec_jx = (int)(sf_d(adx)/dx); rec_ox = (int)((sf_o(adx)-ox)/dx);
        rec_ny = sf_n(ady); rec_jy = (int)(sf_d(ady)/dy); rec_oy = (int)((sf_o(ady)-oy)/dy);

        Fimg = sf_output("out");
        sf_settype(Fimg,SF_COMPLEX);
        sf_oaxa(Fimg,az,1);
        sf_oaxa(Fimg,ax,2);
        sf_oaxa(Fimg,ay,3);
    } else {
        Fdat = sf_output("out");
        par.get("rec_ox",rec_ox,0);
        par.get("rec_oy",rec_oy,0);
        par.get("rec_jt",rec_jt,1);
        par.get("rec_jx",rec_jx,1);
        par.get("rec_jy",rec_jy,1);
        par.get("rec_nt",rec_nt,(nt-1)/rec_jt+1);
        par.get("rec_nx",rec_nx,nx);
        par.get("rec_ny",rec_ny,nx);

        sf_warning("Output data dimensions...");
        adt = sf_maxa(rec_nt,0.,rec_jt*dt); sf_setlabel(adt,"t");        if(verb) sf_raxa(adt);
        adx = sf_maxa(rec_nx,rec_ox*dx,rec_jx*dx); sf_setlabel(adx,"x"); if(verb) sf_raxa(adx);
        ady = sf_maxa(rec_ny,rec_oy*dy,rec_jy*dy); sf_setlabel(ady,"y"); if(verb) sf_raxa(ady);
        sf_settype(Fdat,SF_COMPLEX);
        sf_oaxa(Fdat,adt,1);
        sf_oaxa(Fdat,adx,2);
        sf_oaxa(Fdat,ady,3);
    }
    // now change to padded grid
    rec_dep+= nb;
    rec_ox += nb;
    rec_oy += nb;

    /* create data (and image) array */
    sf_complex ***dat = sf_complexalloc3(rec_nt,rec_nx,rec_ny);
    sf_complex ***img = NULL;
    if (migr) {
        sf_complexread(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);
        img = sf_complexalloc3(nz,nx,ny);
    }

    // source location
    int sou_x; par.get("sou_x",sou_x); // source position x
    int sou_y; par.get("sou_y",sou_y); // source position y
    int sou_z; par.get("sou_z",sou_z); // source position z
    // now change to padded grid
    sou_x += nb;
    sou_y += nb;
    sou_z += nb;
    // mutting first arrival
    bool mute;
    float sou_t0, vel_w;
    if(!migr) {
        par.get("mute",mute,false);
        if(mute) {
            par.get("sou_t0",sou_t0,0.);  // source delay
            par.get("vel_w",vel_w,1500.); // water velocity
        }
    }

    // wavefield snapshot
    bool snap; par.get("snap",snap,false);
    int jsnap; par.get("jsnap",jsnap,nt);

    rtm3d rtm = rtm3d_init(sou_x, sou_y, sou_z, nt, dt, sou_t0, vel_w, roll, rec_dep, rec_ox, rec_oy, rec_jt, rec_jx, rec_jy, rec_nt, rec_nx, rec_ny, snap, jsnap);

    // Gaussian bell
    int nbell; par.get("nbell",nbell,1); // source position z
    if(nbell) bel3d_init(nbell,fdm,rtm);

    sf_complex ***u = sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz)                   \
            shared(u,fdm)
#endif
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
                u[iy][ix][iz] = sf_cmplx(0,0);
            }
        }
    }
    sf_complex ***bu = NULL;
    if(migr) {
        bu = sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz)                   \
            shared(bu,fdm)
#endif
        for        (iy=0; iy<fdm->nypad; iy++) {
            for    (ix=0; ix<fdm->nxpad; ix++) {
                for(iz=0; iz<fdm->nzpad; iz++) {
                    bu[iy][ix][iz] = sf_cmplx(0,0);
                }
            }
        }
    }

    sf_complex ***uc = NULL;
    sf_axis acz = NULL, acx = NULL, acy = NULL;
    sf_file Fwfl = NULL; /* wavefield */
    if(snap) {
        Fwfl = sf_output("wfl"); /* wavefield */
        int   nqz; par.get("nqz",nqz,nz);
        int   nqx; par.get("nqx",nqx,nx);
        int   nqy; par.get("nqy",nqy,ny);
        float oqz; par.get("oqz",oqz,oz);
        float oqx; par.get("oqx",oqx,ox);
        float oqy; par.get("oqy",oqy,oy);

        uc = sf_complexalloc3(nqz,nqx,nqy);

        sf_warning("Output snapshot dimensions...");
        acz = sf_maxa(nqz,oqz,dz); sf_setlabel(acz,"z"); if(verb) sf_raxa(acz);
        acx = sf_maxa(nqx,oqx,dx); sf_setlabel(acx,"x"); if(verb) sf_raxa(acx);
        acy = sf_maxa(nqy,oqy,dy); sf_setlabel(acy,"y"); if(verb) sf_raxa(acy);
        sf_setn(at,rtm->snp_nt); sf_setd(at,rtm->snp_dt); if(verb) sf_raxa(at);

        sf_settype(Fwfl,SF_COMPLEX);
        sf_oaxa(Fwfl,acz,1);
        sf_oaxa(Fwfl,acx,2);
        sf_oaxa(Fwfl,acy,3);
        sf_oaxa(Fwfl,at, 4);
    }

    /*************************************************************/
    /* lowrank onestep propagator */
    lrk3d lrk = lrk3d_init(nr, lt, rt, fdm, dft);

    /*************************************************************/
    /* revolve (checkpointing) parameters */
    if (migr) { /* migration with optimal checkpointing */

        int check, capo, fine, steps, snaps, oldcapo;
        int info;
        enum action whatodo;
        capo = 0;
        snaps = 0;
        steps = nt;
        //par.get("steps",steps); /* maximum num of steps without saving */
        //snaps = adjust(steps);
        par.get("snaps",snaps); /* maximum num of snapshots allowed to be saved */
        par.get("info",info);   /* verbosity of output info about revolve */
        sf_warning("RTM with checkpointing: steps=%d, snaps=%d, info=%d",steps,snaps,info);
        fine = steps + capo;
        check = -1;             /* Neccessary for first call */

        sf_complex **ustor = sf_complexalloc2(nm,snaps); /* allocation memory */
        do {
            oldcapo = capo;
            whatodo = revolve(&check, &capo, &fine, snaps, &info);
            switch(whatodo) {
                case takeshot:
                    if(info > 1) sf_warning(" takeshot at %6d ",capo);
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz,i)                 \
            shared(u,ustor,fdm,check)
#endif
                    for        (iy=0; iy<fdm->nypad; iy++) {
                        for    (ix=0; ix<fdm->nxpad; ix++) {
                            for(iz=0; iz<fdm->nzpad; iz++) {
                                i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz;
                                ustor[check][i] = u[iy][ix][iz];
                            }
                        }
                    }
                    break;
                case advance:
                    if(info > 2) sf_warning(" advance to %7d (prop from %d to %d) ",capo,oldcapo,capo);
                    for(it=oldcapo; it<capo; it++) {
                        /* inject source */
                        inject_bell_src(u, ww[it], rtm);
                        /* forward prop */
                        forward(u, fdm, dft, lrk, spo);
                    }
                    break;
                case firsturn:
                    if(info > 2) sf_warning(" firsturn at %6d ",capo);
                    /* initialize image */
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz)                   \
            shared(img,ny,nx,nz)
#endif
                    for        (iy=0; iy<ny; iy++) {
                        for    (ix=0; ix<nx; ix++) {
                            for(iz=0; iz<nz; iz++) {
                                img[iy][ix][iz] = sf_cmplx(0,0);
                            }
                        }
                    }
                    /* inject data */
                    inject3d(bu, dat, capo, rtm);
                    /* backward prop */
                    reverse(bu, fdm, dft, lrk, spo);
                    /* cross-correlation imaging condition */
                    ccr(img, u, bu, fdm);
                    break;
                case youturn:
                    if(info > 2) sf_warning(" youturn at %7d ",capo);
                    /* inject data */
                    inject3d(bu, dat, capo, rtm);
                    /* backward prop */
                    reverse(bu, fdm, dft, lrk, spo);
                    /* cross-correlation imaging condition */
                    ccr(img, u, bu, fdm);
                    break;
                case restore:
                    if(info > 2) sf_warning(" restore at %7d ",capo);
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz,i)                 \
            shared(u,ustor,fdm,check)
#endif
                    for        (iy=0; iy<fdm->nypad; iy++) {
                        for    (ix=0; ix<fdm->nxpad; ix++) {
                            for(iz=0; iz<fdm->nzpad; iz++) {
                                i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz;
                                u[iy][ix][iz] = ustor[check][i];
                            }
                        }
                    }
                    break;
                case error:
                    sf_warning(" irregular termination of revolve ");
                    switch(info) {
                        case 10: sf_warning(" number of checkpoints stored exceeds checkup, ");
                                 sf_warning(" increase constant 'checkup' and recompile ");
                                 break;
                        case 11: sf_warning(" number of checkpoints stored = %d exceeds snaps = %d, ",check+1,snaps);
                                 sf_warning(" ensure 'snaps' > 0 and increase initial 'fine' ");
                                 break;
                        case 12: sf_warning(" error occurs in numforw ");
                                 break;
                        case 13: sf_warning(" enhancement of 'fine', 'snaps' checkpoints stored, ");
                                 sf_warning(" increase 'snaps' ");
                                 break;
                        case 14: sf_warning(" number of snaps exceeds snapsup, ");
                                 sf_warning(" increase constant 'snapsup' and recompile ");
                                 break;
                        case 15: sf_warning(" number of reps exceeds repsup, ");
                                 sf_warning(" increase constant 'repsup' and recompile ");       
                    }
                    sf_error(" exiting...");
                    break;
                default:
                    break;
            }
        } while((whatodo != terminate) && (whatodo != error));

        /* output image */
        sf_complexwrite(img[0][0],nz*nx*ny,Fimg);

    } else { /* forward modeling of shot records */

        for (it=0; it<nt; it++) {
            if(verb) sf_warning("it=%d/%d;", it, nt-1);

            /* take snapshot */
            if(snap && it%jsnap==0) {
                cut3d_complex(u,uc,fdm,acz,acx,acy);
                sf_complexwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
            }

            /* record data */
            if(it%rec_jt==0) {
                extract3d(u, dat, it, rtm);
            }

            /* inject source */
            inject_bell_src(u, ww[it], rtm);

            /* forward prop */
            forward(u, fdm, dft, lrk, spo);

        }
        if(verb) sf_warning(".");

        /* mute first arrival */
        if(mute) mute3d(dat, fdm, rtm);

        /* output data */
        sf_complexwrite(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);

    }

    /*************************************************************/
    /* clean up memory variables */

    dft3d_finalize();
    lrk3d_finalize();
    if(nbell) bel3d_finalize();

    free(fdm);
    free(dft);
    free(rtm);
    free(lrk);

    free(ww);
    free(**tt); free(*tt); free(tt);
    free(**vel); free(*vel); free(vel);
    free(*lt); free(lt);
    free(*rt); free(rt);
    free(**dat); free(*dat); free(dat);
    free(**u); free(*u); free(u);
    if(migr) {
        free(**img); free(*img); free(img);
        free(**bu); free(*bu); free(bu);
    }
    if(snap) {
        free(**uc); free(*uc); free(uc);
    }

    return 0;
}
