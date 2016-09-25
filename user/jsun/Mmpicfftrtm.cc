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
#include <mpi.h>
#include "lowrank.hh"
/* head files aumatically produced from C programs */
extern "C"{
#include "rtmutil.h"
#include "revolve.h"
}

//using namespace std;

#define GEO 14

int main(int argc, char** argv)
{   
    /*MPI*/
    int cpuid, numprocs, provided;
    //sf_complex *sendbuf, *recvbuf;
    float *sendbuf, *recvbuf;

    //MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided<MPI_THREAD_FUNNELED) sf_warning("MPI does not support OpenMP");

    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    sf_init(argc,argv); // Initialize RSF

    /*************************************************************/
    /* GLOBAL DEFINITION (valid thru entire shot-loop)           */
    /*************************************************************/

    iRSF par(0);

    /* controlling flags */
    bool verb; par.get("verb",verb,false); // verbosity
    bool migr; par.get("migr",migr,false); // adjoint(migration) flag
    bool roll; par.get("roll",roll,false); // rolling v.s. fixed-spread acquisition geometry
    bool dabc; par.get("dabc",dabc,false); // absorbing boundary
    bool snap; par.get("snap",snap,false); // output wavefield snapshots
    bool mute; par.get("mute",mute,false); // mute first arrival (modeling or imaging)
    bool sill; par.get("sill",sill,false); // source illumination for rtm

    /* dimension related */
    int nb;    par.get("nb",nb,0);          // abc width
    if(nb==0) dabc=false;
    float cb;    par.get("cb",cb,1.0);        // abc strength
    int nbell; par.get("nbell",nbell,1);  // source position z
    int jsnap; par.get("jsnap",jsnap,100); // snapshot interval

    /* lowrank parameters */
    int jump;  par.get("jump",jump,1);          // subsampling rate for lowrank decomposition
    int seed;  par.get("seed",seed,time(NULL)); // seed for random number generator
    int npk;   par.get("npk",npk,20);           // maximum rank
    float eps; par.get("eps",eps,1.e-4);        // tolerance/accuracy

    /* media */
    int   media; par.get("media",media,0); // media: 0-> iso, 1-> tti
    int   taper; par.get("taper",taper,0); // tapering interval for tti
    float thres; par.get("thres",thres,1); // tapering threshold for tti
    if (media==1) {
        if (taper==0) taper=4;
        if (thres==1) thres=0.92;
    }

    /* mutting parameters */
    float sou_t0; par.get("sou_t0",sou_t0,0.);  // source delay
    float vel_w;  par.get("vel_w",vel_w,1500.); // water velocity

    /* check pointing parameters */
    int revolve_snaps; par.get("revolve_snaps",revolve_snaps,64); // maximum num of snapshots allowed to be saved
    int info;          par.get("info",info,0);                    // verbosity of output info about revolve

    /* data directory */
    char *dat_dir = sf_getstring("dat_dir");
    if(NULL==dat_dir) sf_error("Need dat_dir=");
    char *wfl_dir = sf_getstring("wfl_dir");
    if(snap && NULL==wfl_dir) sf_error("Need wfl_dir=");
    char *img_dir;
    if (migr) {
    img_dir = sf_getstring("img_dir");
    if(NULL==img_dir) sf_error("Need img_dir=");
    }

    /* file name for geo par */
    sf_file geo = sf_input("geo");  // geometry: sou_z, sou_x, sou_y, rec_dep, rec_ox, rec_oy, rec_nx, rec_ny;
    int esize,sht_num_total,ng1;
    if(SF_INT != sf_gettype(geo)) sf_error("Need int in geo!");
    if(!sf_histint(geo,"esize",&esize)) esize = sizeof("int");
    if(!sf_histint(geo,"n1",&ng1)) sf_error("No n1= in geo!");
    if(ng1!=GEO) sf_error("n1 in geo should be %d",GEO);
    if(!sf_histint(geo,"n2",&sht_num_total)) sf_error("No n2= in geo!");
    if(cpuid==0) sf_warning("A total of %d shots!",sht_num_total);

    /* shot control parameter */
    int sht_set; par.get("sht_set",sht_set,0);            // starting shot index
    int sht_num; par.get("sht_num",sht_num,sht_num_total);// shot number to process

    /* global spatial and time axis */
    sf_file Fwav = sf_input("wav"); // source wavelet
    sf_file Fvel=NULL, Fvelz=NULL, Feta=NULL, Ftheta=NULL; 
    switch(media) {
        case 1:
            Fvel   = sf_input("velx");  // horizontal velocity
            Fvelz  = sf_input("velz");  // vertical velocity
            Feta   = sf_input("eta");   // anelliptic parameter
            Ftheta = sf_input("theta"); // tilt angle (in degrees)
            break;
        default:
            Fvel = sf_input("vel"); // velocity
            break;
    }

    if(verb) sf_warning("Global model dimensions...");
    sf_axis az_g = sf_iaxa(Fvel,1); sf_setlabel(az_g,"z"); if(verb) sf_raxa(az_g); /* depth z */
    sf_axis ax_g = sf_iaxa(Fvel,2); sf_setlabel(ax_g,"x"); if(verb) sf_raxa(ax_g); /* space x */
    sf_axis ay_g = sf_iaxa(Fvel,3); sf_setlabel(ay_g,"y"); if(verb) sf_raxa(ay_g); /* space y */
    sf_axis at_g = sf_iaxa(Fwav,1); sf_setlabel(at_g,"t"); if(verb) sf_raxa(at_g); /* time  t */

    int nz_g = sf_n(az_g); float dz_g = sf_d(az_g); float oz_g = sf_o(az_g);
    int nx_g = sf_n(ax_g); float dx_g = sf_d(ax_g); float ox_g = sf_o(ax_g);
    int ny_g = sf_n(ay_g); float dy_g = sf_d(ay_g); float oy_g = sf_o(ay_g);
    int nt_g = sf_n(at_g); float dt_g = sf_d(at_g); /*float ot_g = sf_o(at_g);*/

    if(cpuid==0) {
        if(ny_g>1) sf_warning("3D Modeling/Migration");
        else sf_warning("2D Modeling/Migration");
    }

    fdm3d fdm_g = fdutil3d_init(verb,false,az_g,ax_g,ay_g,0,1);

    /* source wavelet */
    sf_complex *ww = sf_complexalloc(nt_g); // fast axis: n_time > n_shot
    sf_complexread(ww,nt_g,Fwav);

    /* model parameter */
    float ***vel_g=NULL, ***velx_g=NULL, ***velz_g=NULL, ***eta_g=NULL, ***theta_g=NULL;
    switch(media) {
        case 1:
            velx_g  = sf_floatalloc3(nz_g,nx_g,ny_g); 
            sf_floatread( velx_g[0][0],nz_g*nx_g*ny_g,Fvel  );

            velz_g  = sf_floatalloc3(nz_g,nx_g,ny_g); 
            sf_floatread( velz_g[0][0],nz_g*nx_g*ny_g,Fvelz );

            eta_g   = sf_floatalloc3(nz_g,nx_g,ny_g); 
            sf_floatread(  eta_g[0][0],nz_g*nx_g*ny_g,Feta  );

            theta_g = sf_floatalloc3(nz_g,nx_g,ny_g); 
            sf_floatread(theta_g[0][0],nz_g*nx_g*ny_g,Ftheta);
            break;
        default:
            vel_g = sf_floatalloc3(nz_g,nx_g,ny_g); 
            sf_floatread(vel_g[0][0],nz_g*nx_g*ny_g,Fvel);
            break;
    }

    /* loop index */
    int is,it,iy,ix,iz,i;

    /* output image */
    sf_file Fimage = NULL;
    //sf_complex ***img_g = NULL;
    float ***img_g = NULL;
    if(migr) {
        //img_g = sf_complexalloc3(nz_g,nx_g,ny_g);
        //setval_complex(img_g[0][0],nz_g*nx_g*ny_g,sf_cmplx(0,0));
        img_g = sf_floatalloc3(nz_g,nx_g,ny_g);
        setval(img_g[0][0],nz_g*nx_g*ny_g,0.f);
        if(cpuid==0) {
            Fimage = sf_output("image");
            //sf_settype(Fimage,SF_COMPLEX);
            sf_settype(Fimage,SF_FLOAT);
            sf_oaxa(Fimage,az_g,1);
            sf_oaxa(Fimage,ax_g,2);
            sf_oaxa(Fimage,ay_g,3);
        }
    }

    /*************************************************************/
    /* LOCAL DEFINITION (shot/receiver number dependent)         */
    /*************************************************************/

    /* determine shot number */
    for (is=0; is*numprocs<sht_num; is++) {
        int sht_id = is*numprocs + cpuid + sht_set;
        if (sht_id < sht_set+sht_num) { // do the work, otherwise idle
            sf_warning("Node #%d working on shot #%d",cpuid,sht_id);

            /* file name for shot record */
            char sht_id_c[10];
            itoa(sht_id,sht_id_c);
            char dat_file[100];
            strcpy(dat_file,dat_dir);
            strcat(dat_file,"/shot-");
            strcat(dat_file,sht_id_c);
            strcat(dat_file,".rsf");
            if(verb) sf_warning("cpuid=%d, data file=%s",cpuid,dat_file);

            /* read geometry file and prepare data */
            int mod_oz_g, mod_ox_g, mod_oy_g, mod_nz_g, mod_nx_g, mod_ny_g, sou_z_g, sou_x_g, sou_y_g, rec_dep_g, rec_ox_g, rec_oy_g, rec_jt_g, rec_jx_g, rec_jy_g, rec_nt_g, rec_nx_g, rec_ny_g;
            sf_seek(geo,sht_id*esize*GEO,SEEK_SET); // geometry: mod_oz, mod_ox, mod_oy, mod_nz, mod_nx, mod_ny, sou_z, sou_x, sou_y, rec_dep, rec_ox, rec_oy, rec_nx, rec_ny;
            sf_intread(&mod_oz_g, 1,geo);
            sf_intread(&mod_ox_g, 1,geo);
            sf_intread(&mod_oy_g, 1,geo);
            sf_intread(&mod_nz_g, 1,geo);
            sf_intread(&mod_nx_g, 1,geo);
            sf_intread(&mod_ny_g, 1,geo);
            sf_intread(&sou_z_g,  1,geo);
            sf_intread(&sou_x_g,  1,geo);
            sf_intread(&sou_y_g,  1,geo);
            sf_intread(&rec_dep_g,1,geo);
            sf_intread(&rec_ox_g, 1,geo);
            sf_intread(&rec_oy_g, 1,geo);
            sf_intread(&rec_nx_g, 1,geo);
            sf_intread(&rec_ny_g, 1,geo);

            rec_jt_g = 1; rec_jx_g = 1; rec_jy_g = 1; rec_nt_g = nt_g; // these pars can be read from input as well

            if(ny_g>1) {
                if( rec_ox_g<0 || rec_ox_g+rec_nx_g*rec_jx_g>nx_g || rec_oy_g<0 || rec_oy_g+rec_ny_g*rec_jy_g>ny_g ) sf_error("Data dimension out of bound!");
            } else {
                if( rec_ox_g<0 || rec_ox_g+rec_nx_g*rec_jx_g>nx_g ) sf_error("Data dimension out of bound!");
                if( mod_oy_g!=0)  sf_error("Set mod_oy_g to ZERO for 2D jobs!");
                if( mod_ny_g!=1)  sf_error("Set mod_ny_g to ONE for 2D jobs!");
                if( sou_y_g!=0)  sf_error("Set sou_y_g to ZERO for 2D jobs!");
                if( rec_oy_g!=0) sf_error("Set rec_oy_g to ZERO for 2D jobs!");
                if( rec_ny_g!=1) sf_error("Set rec_ny_g to ONE for 2D jobs!");
            }

            sf_file Fdat=NULL;
            sf_axis adt = NULL, adx = NULL, ady = NULL;
            if (migr) {
                Fdat = sf_input(dat_file);
                if(verb) sf_warning("Input data dimensions...");
                adt = sf_iaxa(Fdat,1); if(verb) sf_raxa(adt); /* time  t */
                adx = sf_iaxa(Fdat,2); if(verb) sf_raxa(adx); /* space x */
                ady = sf_iaxa(Fdat,3); if(verb) sf_raxa(ady); /* space y */
                if (rec_nt_g != sf_n(adt) || rec_jt_g != (int)(sf_d(adt)/dt_g+0.5)) /*rec_ot_g = (int)((sf_o(adt)-ot_g)/dt_g);*/
                    sf_error("Input data first dimension unmatched! %d, %d, %d, %d",rec_nt_g,sf_n(adt),rec_jt_g,(int)(sf_d(adt)/dt_g+0.5));
                if (rec_nx_g != sf_n(adx) || rec_jx_g != (int)(sf_d(adx)/dx_g+0.5) || rec_ox_g != (int)((sf_o(adx)-ox_g)/dx_g+0.5))
                    sf_error("Input data second dimension unmatched! %d, %d, %d, %d, %d, %d", rec_nx_g,sf_n(adx),rec_jx_g,(int)(sf_d(adx)/dx_g+0.5),rec_ox_g,(int)((sf_o(adx)-ox_g)/dx_g+0.5));
                if (rec_ny_g != sf_n(ady) || rec_jy_g != (int)(sf_d(ady)/dy_g+0.5) || rec_oy_g != (int)((sf_o(ady)-oy_g)/dy_g+0.5))
                    sf_error("Input data third dimension unmatched! %d, %d, %d, %d, %d, %d",rec_ny_g,sf_n(ady),rec_jy_g,(int)(sf_d(ady)/dy_g+0.5),rec_oy_g,(int)((sf_o(ady)-oy_g)/dy_g+0.5));
            } else {
                Fdat = sf_output(dat_file);
                if(verb) sf_warning("Output data dimensions...");
                adt = sf_maxa(rec_nt_g,0.,rec_jt_g*dt_g);                 sf_setlabel(adt,"t"); if(verb) sf_raxa(adt);
                adx = sf_maxa(rec_nx_g,ox_g+rec_ox_g*dx_g,rec_jx_g*dx_g); sf_setlabel(adx,"x"); if(verb) sf_raxa(adx);
                ady = sf_maxa(rec_ny_g,oy_g+rec_oy_g*dy_g,rec_jy_g*dy_g); sf_setlabel(ady,"y"); if(verb) sf_raxa(ady);
                //sf_settype(Fdat,SF_COMPLEX);
                sf_settype(Fdat,SF_FLOAT);
                sf_oaxa(Fdat,adt,1);
                sf_oaxa(Fdat,adx,2);
                sf_oaxa(Fdat,ady,3);
            }

            /* calculate local coordinate */
            int nt = nt_g;     float dt = dt_g;
            int nz = mod_nz_g; float dz = dz_g; float oz = oz_g + mod_oz_g*dz_g;
            int nx = mod_nx_g; float dx = dx_g; float ox = ox_g + mod_ox_g*dx_g;
            int ny = mod_ny_g; float dy = dy_g; float oy = oy_g + mod_oy_g*dy_g;

            int sou_z = sou_z_g-mod_oz_g;     int sou_x = sou_x_g-mod_ox_g;   int sou_y = sou_y_g-mod_oy_g;
            int rec_dep = rec_dep_g-mod_oz_g; int rec_ox = rec_ox_g-mod_ox_g; int rec_oy = rec_oy_g-mod_oy_g;
            int rec_jt = rec_jt_g;            int rec_jx = rec_jx_g;          int rec_jy = rec_jy_g;
            int rec_nt = rec_nt_g;            int rec_nx = rec_nx_g;          int rec_ny = rec_ny_g;

            if(verb) sf_warning("Local model dimension...");
            sf_axis az = sf_maxa(nz,oz,dz); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
            sf_axis ax = sf_maxa(nx,ox,dx); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);
            sf_axis ay = sf_maxa(ny,oy,dy); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay);

            /* change to padded grid */
            sou_z += nb; sou_x += nb; sou_y += (ny>1)?nb:0;
            rec_dep+= nb; rec_ox += nb; rec_oy += (ny>1)?nb:0;

            /* space domain prep */
            fdm3d fdm = fdutil3d_init(verb,false,az,ax,ay,nb,1);
            int nm = fdm->nzpad*fdm->nxpad*fdm->nypad;

            float ***vel=NULL, ***velx=NULL, ***velz=NULL, ***eta=NULL, ***theta=NULL;
            float ***tmp = sf_floatalloc3(nz,nx,ny);
            switch(media) {
                case 1:
                    velx = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                    cut3d(velx_g,tmp,fdm_g,az,ax,ay);
                    if(ny>1) expand3d(tmp,velx,fdm);
                    else expand2d(tmp[0],velx[0],fdm);

                    velz = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                    cut3d(velz_g,tmp,fdm_g,az,ax,ay);
                    if(ny>1) expand3d(tmp,velz,fdm);
                    else expand2d(tmp[0],velz[0],fdm);


                    eta = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                    cut3d(eta_g,tmp,fdm_g,az,ax,ay);
                    if(ny>1) expand3d(tmp,eta,fdm);
                    else expand2d(tmp[0],eta[0],fdm);

                    theta = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                    cut3d(theta_g,tmp,fdm_g,az,ax,ay);
                    if(ny>1) expand3d(tmp,theta,fdm);
                    else expand2d(tmp[0],theta[0],fdm);
                    break;
                default:
                    vel = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                    cut3d(vel_g,tmp,fdm_g,az,ax,ay);
                    if(ny>1) expand3d(tmp,vel,fdm);
                    else expand2d(tmp[0],vel[0],fdm);
                    break;
            }

            /* set up abc */
            sponge spo=NULL;
            if (dabc) {
                spo = sponge_make2(fdm->nb,cb);
            }

            /* wavenumber domain prep */
            dft3d dft = dft3d_init(1,false,false,fdm);
            int nn = dft->nkz*dft->nkx*dft->nky;

            /* set up tapering */
            if (taper) tap3d_init(thres,dft);

            /*************************************************************/
            /* perform lowrank decomposition on-the-fly */
            sf_warning("Lowrank decomposition...");
            lowrank_init(jump, seed, npk, eps, dt, media, fdm, dft);
            switch(media) {
                case 1:
                    lowrank_tti(velx[0][0],velz[0][0],eta[0][0],theta[0][0]);
                    break;
                default:
                    lowrank_iso(vel[0][0]);
                    break;
            }
            int nr = lowrank_rank();
            sf_complex **lt = sf_complexalloc2(nm,nr);
            sf_complex **rt = sf_complexalloc2(nr,nn);
            lowrank_mat(lt,rt);

            /* create data (and image) array */
            //sf_complex ***dat = sf_complexalloc3(rec_nt,rec_nx,rec_ny);
            //sf_complex ***img = NULL, ***sil = NULL;
            float ***dat = sf_floatalloc3(rec_nt,rec_nx,rec_ny);
            float ***img = NULL, ***sil = NULL;

            rtm3d rtm = rtm3d_init(sou_x, sou_y, sou_z, nt, dt, sou_t0, vel_w, roll, rec_dep, rec_ox, rec_oy, rec_jt, rec_jx, rec_jy, rec_nt, rec_nx, rec_ny, snap, jsnap);

            // Gaussian bell
            if(nbell) bel3d_init(nbell,fdm,rtm);

            // state and adjoint wavefield
            sf_complex ***u = sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
            setval_complex(u[0][0],fdm->nzpad*fdm->nxpad*fdm->nypad,sf_cmplx(0,0));
            sf_complex ***bu = NULL;

            // output file
            int   nqz=nz; int   nqx=nx; int   nqy=ny;
            float oqz=oz; float oqx=ox; float oqy=oy;

            if(verb) sf_warning("Output shot-img dimensions...");
            sf_axis acz = sf_maxa(nqz,oqz,dz); sf_setlabel(acz,"z"); if(verb) sf_raxa(acz);
            sf_axis acx = sf_maxa(nqx,oqx,dx); sf_setlabel(acx,"x"); if(verb) sf_raxa(acx);
            sf_axis acy = sf_maxa(nqy,oqy,dy); sf_setlabel(acy,"y"); if(verb) sf_raxa(acy);

            sf_file Fimg = NULL; /* single shot image */
            if(migr) {
                char img_file[100];
                strcpy(img_file,img_dir);
                strcat(img_file,"/img-");
                strcat(img_file,sht_id_c);
                strcat(img_file,".rsf");
                if(verb) sf_warning("cpuid=%d, single shot image file=%s",cpuid,img_file);

                Fimg = sf_output(img_file);
                //sf_settype(Fimg,SF_COMPLEX);
                sf_settype(Fimg,SF_FLOAT);
                sf_oaxa(Fimg,acz,1);
                sf_oaxa(Fimg,acx,2);
                sf_oaxa(Fimg,acy,3);
            }

            sf_complex ***uc = NULL;
            sf_axis act = NULL;
            sf_file Fwfl = NULL; /* wavefield */
            if(snap) {
                uc = sf_complexalloc3(nqz,nqx,nqy);

                if(verb) sf_warning("Output snapshot dimensions...");
                act = sf_maxa(rtm->snp_nt,0,rtm->snp_dt); sf_setlabel(act,"t"); if(verb) sf_raxa(act);

                char wfl_file[100];
                strcpy(wfl_file,wfl_dir);
                strcat(wfl_file,"/wfl-");
                strcat(wfl_file,sht_id_c);
                strcat(wfl_file,".rsf");
                if(verb) sf_warning("cpuid=%d, wavefield file=%s",cpuid,wfl_file);

                Fwfl = sf_output(wfl_file);
                sf_settype(Fwfl,SF_COMPLEX);
                sf_oaxa(Fwfl,acz,1);
                sf_oaxa(Fwfl,acx,2);
                sf_oaxa(Fwfl,acy,3);
                sf_oaxa(Fwfl,act,4);
            }

            /*************************************************************/
            /* lowrank onestep propagator */
            lrk3d lrk = lrk3d_init(nr, lt, rt, fdm, dft);

            /*************************************************************/
            /* revolve (checkpointing) parameters */
            if (migr) { /* migration with optimal checkpointing */

                int check, capo, fine, steps, oldcapo;
                enum action whatodo;
                capo = 0;
                steps = nt;
                //revolve_snaps = adjust(steps);
                revolve_snaps = (revolve_snaps>nt)? nt:revolve_snaps;
                if(verb) sf_warning("RTM with checkpointing: steps=%d, snaps=%d, info=%d",steps,revolve_snaps,info);
                fine = steps + capo;
                check = -1;             /* Neccessary for first call */

                sf_complex **ustor = sf_complexalloc2(nm,revolve_snaps); /* allocation memory */
                do {
                    oldcapo = capo;
                    whatodo = revolve(&check, &capo, &fine, revolve_snaps, &info);
                    switch(whatodo) {
                        case takeshot:
                            if(info > 1) sf_warning("node#%3d takeshot at %6d ",cpuid,capo);
#ifdef _OPENMP
#pragma omp parallel for       \
                            private(iy,ix,iz,i)                \
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
                            if(info > 2) sf_warning("node#%3d advance to %7d (prop from %d to %d) ",cpuid,capo,oldcapo,capo);
                            for(it=oldcapo; it<capo; it++) {
                                /* inject source */
                                inject_bell_src(u, ww[it], rtm);
                                /* forward prop */
                                forward(u, taper!=0 && it%taper==0, fdm, dft, lrk, spo);
                            }
                            break;
                        case firsturn:
                            if(info > 2) sf_warning("node#%3d firsturn at %6d ",cpuid,capo);
                            /* read data */
                            //sf_complexread(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);
                            sf_floatread(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);
                            /* mute first arrival */
                            if(mute) mute3d(dat, fdm, rtm);
                            /* initialize image */
                            //img = sf_complexalloc3(nz,nx,ny);
                            //setval_complex(img[0][0],nz*nx*ny,sf_cmplx(0,0));
                            img = sf_floatalloc3(nz,nx,ny);
                            setval(img[0][0],nz*nx*ny,0.f);
                            if(sill) {
                                //sil = sf_complexalloc3(nz,nx,ny);
                                //setval_complex(sil[0][0],nz*nx*ny,sf_cmplx(0,0));
                                sil = sf_floatalloc3(nz,nx,ny);
                                setval(sil[0][0],nz*nx*ny,0.f);
                            }
                            /* initialize adjoint wavefield */
                            bu = sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
                            setval_complex(bu[0][0],fdm->nzpad*fdm->nxpad*fdm->nypad,sf_cmplx(0,0));

                            /* 4 - backward prop */
                            reverse(bu, taper!=0 && capo%taper==0, fdm, dft, lrk, spo);
                            /* 3 - inject data */
                            if(capo%rec_jt==0) {
                                inject3d(bu, dat, capo, rtm);
                            }
                            /* 2 - take snapshot */
                            if(snap && capo%jsnap==0) {
                                cut3d_complex(bu,uc,fdm,acz,acx,acy);
                                sf_complexwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
                            }
                            /* 1 - cross-correlation imaging condition */
                            ccrf(img, u, bu, fdm);
                            if(sill) ccrf(sil, u,  u, fdm);
                            break;
                        case youturn:
                            if(info > 2) sf_warning("node#%3d youturn at %7d ",cpuid,capo);
                            /* 4 - backward prop */
                            reverse(bu, taper!=0 && capo%taper==0, fdm, dft, lrk, spo);
                            /* 3 - inject data */
                            if(capo%rec_jt==0) {
                                inject3d(bu, dat, capo, rtm);
                            }
                            /* 2 - take snapshot */
                            if(snap && capo%jsnap==0) {
                                cut3d_complex(bu,uc,fdm,acz,acx,acy);
                                sf_complexwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
                            }
                            /* 1 - cross-correlation imaging condition */
                            ccrf(img, u, bu, fdm);
                            if(sill) ccrf(sil, u,  u, fdm);
                            break;
                        case restore:
                            if(info > 2) sf_warning("node#%3d restore at %7d ",cpuid,capo);
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
                                case 11: sf_warning(" number of checkpoints stored = %d exceeds snaps = %d, ",check+1,revolve_snaps);
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

                if(sill) {
#ifdef _OPENMP
#pragma omp parallel for                    \
                    private(iy,ix,iz)       \
                    shared(img,sil)
#endif
                    for        (iy=0; iy<ny; iy++) {
                        for    (ix=0; ix<nx; ix++) {
                            for(iz=0; iz<nz; iz++) {
                                img[iy][ix][iz] = img[iy][ix][iz]/(sil[iy][ix][iz]+SF_EPS);
//#ifdef SF_HAS_COMPLEX_H
//                                img[iy][ix][iz] = img[iy][ix][iz]/(sil[iy][ix][iz]+SF_EPS);
//#else
//                                img[iy][ix][iz] = sf_cdiv(img[iy][ix][iz],sf_cadd(sil[iy][ix][iz],sf_cmplx(SF_EPS,0)));
//#endif
                            }
                        }
                    }
                }

                /* output image */
                //sf_complexwrite(img[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fimg);
                sf_floatwrite(img[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fimg);

                /* stack image */
#ifdef _OPENMP
#pragma omp parallel for       \
                private(iy,ix,iz)                  \
                shared(img_g,img,ny,nx,nz,mod_oy_g,mod_ox_g,mod_oz_g)
#endif
                for        (iy=0; iy<ny; iy++) {
                    for    (ix=0; ix<nx; ix++) {
                        for(iz=0; iz<nz; iz++) {
                            img_g[iy+mod_oy_g][ix+mod_ox_g][iz+mod_oz_g] += img[iy][ix][iz];
//#ifdef SF_HAS_COMPLEX_H
//                            img_g[iy+mod_oy_g][ix+mod_ox_g][iz+mod_oz_g] += img[iy][ix][iz];
//#else
//                            img_g[iy+mod_oy_g][ix+mod_ox_g][iz+mod_oz_g] = sf_cadd(img_g[iy+mod_oy_g][ix+mod_ox_g][iz+mod_oz_g],img[iy][ix][iz]);
//#endif
                        }
                    }
                }

                free(*ustor); free(ustor);

            } else { /* forward modeling of shot records */

                for (it=0; it<nt; it++) {
                    if(verb) sf_warning("it=%d/%d;", it, nt-1);

                    /* 1 - inject source */
                    inject_bell_src(u, ww[it], rtm);
                    /* 2 - take snapshot */
                    if(snap && it%jsnap==0) {
                        cut3d_complex(u,uc,fdm,acz,acx,acy);
                        sf_complexwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
                    }
                    /* 3 - record data */
                    if(it%rec_jt==0) {
                        extract3d(u, dat, it, rtm);
                    }
                    /* 4 - forward prop */
                    forward(u, taper!=0 && it%taper==0, fdm, dft, lrk, spo);

                }
                if(verb) sf_warning(".");

                /* mute first arrival */
                if(mute) mute3d(dat, fdm, rtm);

                /* output data */
                //sf_complexwrite(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);
                sf_floatwrite(dat[0][0],rec_nt*rec_nx*rec_ny,Fdat);

            }

            /*************************************************************/
            /* clean up memory variables */
            dft3d_finalize();
            tap3d_finalize();
            lrk3d_finalize();
            if(nbell) bel3d_finalize();

            free(fdm);
            free(dft);
            free(rtm);
            free(lrk);
            if(dabc) free(spo);

            free(**tmp); free(*tmp); free(tmp);
            if(NULL!=vel  ) { free(**vel  ); free(*vel  ); free(vel  ); }
            if(NULL!=velx ) { free(**velx ); free(*velx ); free(velx ); }
            if(NULL!=velz ) { free(**velz ); free(*velz ); free(velz ); }
            if(NULL!=eta  ) { free(**eta  ); free(*eta  ); free(eta  ); }
            if(NULL!=theta) { free(**theta); free(*theta); free(theta); }
            free(*lt); free(lt);
            free(*rt); free(rt);
            free(**dat); free(*dat); free(dat);
            free(**u); free(*u); free(u);
            if(migr) {
                free(**img); free(*img); free(img);
                if(sill) { free(**sil); free(*sil); free(sil); }
                free(**bu); free(*bu); free(bu);
            }
            if(snap) {
                free(**uc); free(*uc); free(uc);
            }

            sf_fileclose(Fdat);
            if(migr) sf_fileclose(Fimg);
            if(snap) sf_fileclose(Fwfl);
        } else {
            sf_warning("Node #%d idling...");
        }
    }

    /* gather image */
    if (migr) {
      if (cpuid==0) {
#if MPI_VERSION >= 2
	//sendbuf = (sf_complex *) MPI_IN_PLACE;
	sendbuf = (float *) MPI_IN_PLACE;
#else /* will fail */
	sendbuf = NULL;
#endif 
	recvbuf = img_g[0][0];
      } else {
	sendbuf = img_g[0][0];
      	recvbuf = NULL;
      }
      //MPI_Reduce(sendbuf, recvbuf, nz_g*nx_g*ny_g, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD); 
      MPI_Reduce(sendbuf, recvbuf, nz_g*nx_g*ny_g, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

      /* output image */
      //if (cpuid==0) sf_complexwrite(img_g[0][0],nz_g*nx_g*ny_g,Fimage);
      if (cpuid==0) sf_floatwrite(img_g[0][0],nz_g*nx_g*ny_g,Fimage);
    }

    free(ww);
    if(NULL!=vel_g  ) { free(**vel_g  ); free(*vel_g  ); free(vel_g  ); }
    if(NULL!=velx_g ) { free(**velx_g ); free(*velx_g ); free(velx_g ); }
    if(NULL!=velz_g ) { free(**velz_g ); free(*velz_g ); free(velz_g ); }
    if(NULL!=eta_g  ) { free(**eta_g  ); free(*eta_g  ); free(eta_g  ); }
    if(NULL!=theta_g) { free(**theta_g); free(*theta_g); free(theta_g); }
    if(migr) { free(**img_g); free(*img_g); free(img_g); }

    MPI_Finalize();
    exit(0);
}
