#include <rsf.h>
#include "rtmutil.h"
#include "fftomp.h"

#include "ksutil.h"
/*^*/

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#ifndef _rtmutil_h

typedef struct lrk3 *lrk3d;
/*^*/

typedef struct rtm3 *rtm3d;
/*^*/

struct lrk3{
    bool adj;
    int seed;
    int npk;
    float eps;
    int nrank;
    sf_complex **lt;
    sf_complex **rt;
};
/*^*/

struct rtm3{
    int sou_x;
    int sou_y;
    int sou_z;
    int sou_nt;
    float sou_dt;
    float sou_t0;
    float vel_w;
    bool roll;
    int rec_dep;
    int rec_ox;
    int rec_oy;
    int rec_nt;
    int rec_nx;
    int rec_ny;
    int rec_jt;
    int rec_jx;
    int rec_jy;
    int sht_num_pad;
    int sht_num_ori;
    int sht_id_cur;
    bool snap;
    int snp_nt;
    int snp_jt;
    float snp_dt;
};
/*^*/

#endif

/*------------------------------------------------------------*/
static sf_complex *cwave =NULL;
static sf_complex *cwavem=NULL;
static sf_complex *wavem =NULL;
static sf_complex **waves=NULL;

lrk3d lrk3d_init(int n2,
                 sf_complex **lft,
                 sf_complex **rht,
                 fdm3d fdm,
                 dft3d dft)
/*< prepare lowrank arrays for rite method >*/
{
    lrk3d lrk;
    int nk,nxyz;

    lrk = (lrk3d) sf_alloc(1,sizeof(*lrk));
    lrk->nrank=n2;
    lrk->lt=lft;
    lrk->rt=rht;

    nxyz = fdm->nypad*fdm->nxpad*fdm->nzpad;
    nk   = dft->nky*dft->nkx*dft->nkz;

    cwave  = sf_complexalloc(nk);
    wavem  = sf_complexalloc(nxyz);
    //cwaves = sf_complexalloc2(nk,lrk->nrank);
    cwavem = sf_complexalloc(nk);
    //waves  = sf_complexalloc2(nxyz,lrk->nrank);
    waves  = sf_complexalloc2(nk,lrk->nrank); /* double usage :) */

    return lrk;
}

/*------------------------------------------------------------*/
void lrk3d_apply(sf_complex *uo,
                 sf_complex *ui,
                 bool adj,
                 fdm3d fdm,
                 dft3d dft,
                 lrk3d lrk)
/*< apply lowrank matrices for time stepping (can be in-place) >*/
{
    int nxyz,nk,ik,im,iz,ix,iy,i;
    sf_complex c;

    nxyz = fdm->nypad*fdm->nxpad*fdm->nzpad;
    nk   = dft->nky*dft->nkx*dft->nkz;

    if (adj) { /* backward propagation - NSPS */

        for (im = 0; im < lrk->nrank; im++) {
#ifdef _OPENMP
#pragma omp parallel for			\
            private(iy,ix,iz,i)                     \
            shared(fdm,lrk,wavem,ui)
#endif
            for (iy=0; iy<fdm->nypad; iy++)         {
                for (ix=0; ix<fdm->nxpad; ix++)     {
                    for (iz=0; iz<fdm->nzpad; iz++) {
                        i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz; /* flattened coordinate */
#ifdef SF_HAS_COMPLEX_H
                        wavem[i] = conjf(lrk->lt[im][i])*ui[i];
#else
                        wavem[i] = sf_cmul(conjf(lrk->lt[im][i]),ui[i]);
#endif
                    }
                }
            }
            fft(wavem,waves[im]);
        }

#ifdef _OPENMP
#pragma omp parallel for              \
        private(im,ik,c)                \
        shared(lrk,cwave,waves,nk)
#endif
        for (ik=0; ik<nk; ik++) {
            c = sf_cmplx(0.,0.);
            for (im=0; im<lrk->nrank; im++) {
#ifdef SF_HAS_COMPLEX_H
                c += waves[im][ik]*conjf(lrk->rt[ik][im]);
#else
                c = sf_cadd(c,sf_cmul(waves[im][ik],conjf(lrk->rt[ik][im])));
#endif
            }
            cwave[ik] = c;
        }

        ifft(uo,cwave);

    } else { /* forward propagation - PSPI */

        fft(ui,cwave);

        for (im=0; im<lrk->nrank; im++) {
#ifdef _OPENMP
#pragma omp parallel for              \
        private(ik)                   \
        shared(lrk,cwavem,cwave,nk)
#endif
            for (ik=0; ik<nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                cwavem[ik] = cwave[ik]*lrk->rt[ik][im];
#else
                cwavem[ik] = sf_cmul(cwave[ik],lrk->rt[ik][im]);
#endif
            }
            ifft(waves[im],cwavem);
        }

#ifdef _OPENMP
#pragma omp parallel for			\
        private(iy,ix,iz,i,c,im)                \
        shared(fdm,lrk,waves,uo)
#endif
        for (iy=0; iy<fdm->nypad; iy++)         {
            for (ix=0; ix<fdm->nxpad; ix++)     {
                for (iz=0; iz<fdm->nzpad; iz++) {
                    i = (iy*fdm->nxpad + ix)*fdm->nzpad + iz; /* flattened coordinate */
                    c = sf_cmplx(0.,0.);
                    for (im = 0; im < lrk->nrank; im++) {
#ifdef SF_HAS_COMPLEX_H
                        c += lrk->lt[im][i]*waves[im][i];
#else
                        c += sf_cmul(lrk->lt[im][i],waves[im][i]);
#endif
                    }
                    uo[i] = c;
                }
            }
        }

    }
}

void lrk3d_finalize()
/*< prepare lowrank arrays for rite method >*/
{
    free(cwave);
    free(wavem);
    free(cwavem);
    free(*waves); free(waves);
}


/*------------------------------------------------------------*/
rtm3d rtm3d_init(int sou_x_,
                 int sou_y_,
                 int sou_z_,
                 int sou_nt_,
                 float sou_dt_,
                 float sou_t0_,
                 float vel_w_,
                 bool roll_,
                 int rec_dep_,
                 int rec_ox_,
                 int rec_oy_,
                 int rec_jt_,
                 int rec_jx_,
                 int rec_jy_,
                 int rec_nt_,
                 int rec_nx_,
                 int rec_ny_,
                 bool snap_,
                 int snp_jt_)
/*< initialize rtm >*/
{
    rtm3d rtm;

    rtm = (rtm3d) sf_alloc(1,sizeof(*rtm));

    rtm->sou_x = sou_x_;
    rtm->sou_y = sou_y_;
    rtm->sou_z = sou_z_;
    rtm->sou_nt = sou_nt_;
    rtm->sou_dt = sou_dt_;
    rtm->sou_t0 = sou_t0_;
    rtm->vel_w = vel_w_;
    rtm->roll = roll_;
    rtm->rec_dep = rec_dep_;
    rtm->rec_ox = rec_ox_;
    rtm->rec_oy = rec_oy_;
    rtm->rec_jt = rec_jt_;
    rtm->rec_jx = rec_jx_;
    rtm->rec_jy = rec_jy_;
    rtm->rec_nt = rec_nt_;
    rtm->rec_nx = rec_nx_;
    rtm->rec_ny = rec_ny_;
    rtm->snap = snap_;
    rtm->snp_jt = snp_jt_;
    rtm->snp_dt = snp_jt_*sou_dt_;
    rtm->snp_nt = 0;
    rtm->snp_nt = (snap_)? (sou_nt_-1)/snp_jt_+1 : 0;

    return rtm;
}

/*------------------------------------------------------------*/
static float *ktp=NULL;
static int nktp;

void tap3d_init(float thres,
                dft3d dft)
/*< init tapering array for tti wave propagation >*/
{
    int iy,ix,iz,nktp,ik;
    float ky,kx,kz,ky_trs,kx_trs,kz_trs,ktmp;

    nktp = dft->nky*dft->nkx*dft->nkz;
    ktp = sf_floatalloc(nktp);

    ky_trs = thres*fabs(dft->oky);
    kx_trs = thres*fabs(dft->okx);
    kz_trs = thres*fabs(dft->okz);
#ifdef _OPENMP
#pragma omp parallel for			\
    private(iy,ix,iz,ik,ktmp,ky,kx,kz)          \
    shared(ktp,ky_trs,kx_trs,kz_trs,dft)
#endif
    for (iy=0; iy<dft->nky; iy++)         {
        ky = dft->oky+iy*dft->dky;
        for (ix=0; ix<dft->nkx; ix++)     {
            kx = dft->okx+ix*dft->dkx;
            for (iz=0; iz<dft->nkz; iz++) {
                kz = dft->okz+iz*dft->dkz;
                ik = (iy*dft->nkx + ix)*dft->nkz + iz; /* flattened coordinate */
                ktmp = 1.;
                ktmp *= (fabs(ky)>ky_trs)? powf((fabs(dft->oky)-fabs(ky)+ky_trs)/dft->oky,2) : 1.;
                ktmp *= (fabs(kx)>kx_trs)? powf((fabs(dft->okx)-fabs(kx)+kx_trs)/dft->okx,2) : 1.;
                ktmp *= (fabs(kz)>kz_trs)? powf((fabs(dft->okz)-fabs(kz)+kz_trs)/dft->okz,2) : 1.;
                ktp[ik] = ktmp;
            }
        }
    }
}

void tap3d_apply(sf_complex *u)
/*< apply tapering to wavefield for tti wave propagation >*/
{
    int ik;

#ifdef _OPENMP
#pragma omp parallel for              \
        private(ik)                   \
        shared(u,ktp,nktp)
#endif
    for (ik=0; ik<nktp; ik++) {
#ifdef SF_HAS_COMPLEX_H
        u[ik] = u[ik]*ktp[ik];
#else
        u[ik] = sf_crmul(u[ik],ktp[ik]);
#endif
    }
}

void tap3d_finalize()
/*< destroy tapering array >*/
{
    free(ktp);
    ktp = NULL;
}

/*------------------------------------------------------------*/
static float***bel3d=NULL;
static int nbel,nbel_y;

void bel3d_init(int n,
                fdm3d fdm,
                rtm3d rtm)
/*< init bell taper >*/
{
    int   iz,ix,iy;
    float s;

    if(fdm->nypad>1) {
        nbel = n;
        nbel_y = n;
        if( rtm->sou_z<(fdm->nb+nbel)   || rtm->sou_z>(fdm->nzpad-fdm->nb-nbel) ||
            rtm->sou_x<(fdm->nb+nbel)   || rtm->sou_x>(fdm->nxpad-fdm->nb-nbel) ||
            rtm->sou_y<(fdm->nb+nbel_y) || rtm->sou_y>(fdm->nypad-fdm->nb-nbel_y) )
        sf_error("Bell taper width too big! sou_z=%d, sou_x=%d, sou_y=%d, nb=%d, nbel=%d",rtm->sou_z,rtm->sou_x,rtm->sou_y,fdm->nb,nbel);
    } else {
        nbel = n;
        nbel_y = 0;
        if( rtm->sou_z<(fdm->nb+nbel) || rtm->sou_z>(fdm->nzpad-fdm->nb-nbel) ||
            rtm->sou_x<(fdm->nb+nbel) || rtm->sou_x>(fdm->nxpad-fdm->nb-nbel) )
        sf_error("Bell taper width too big! sou_z=%d, sou_x=%d, nb=%d, nbel=%d",rtm->sou_z,rtm->sou_x,fdm->nb,nbel);
    }

    s = (nbel==0)? 1 : 2.0/(nbel*nbel);

    bel3d=sf_floatalloc3(2*nbel+1,2*nbel+1,2*nbel_y+1);

    for        (iy=-nbel_y;iy<=nbel_y;iy++) {
        for    (ix=-nbel  ;ix<=nbel  ;ix++) {
            for(iz=-nbel  ;iz<=nbel  ;iz++) {
                bel3d[nbel_y+iy][nbel+ix][nbel+iz] = exp(-(iz*iz+ix*ix+iy*iy)*s);
            }
        }    
    }

}

void inject_bell_src(sf_complex ***u,
                     sf_complex ww,
                     rtm3d rtm)
/*< inject source wavelet into 3d cube >*/
{
    int iy, ix, iz;

    /* inject source */
    for        (iy=-nbel_y;iy<=nbel_y;iy++) {
	for    (ix=-nbel  ;ix<=nbel  ;ix++) {
	    for(iz=-nbel  ;iz<=nbel  ;iz++) {
#ifdef SF_HAS_COMPLEX_H
                u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz] += ww*bel3d[nbel_y+iy][nbel+ix][nbel+iz];
#else
                u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz] = sf_cadd(u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz],sf_crmul(ww,bel3d[nbel_y+iy][nbel+ix][nbel+iz]));
#endif
            }
        }
    }

//    for(iy=-1;iy<=1;iy++) {
//        for(ix=-1;ix<=1;ix++) {
//            for(iz=-1;iz<=1;iz++) {
//#ifdef SF_HAS_COMPLEX_H
//                u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz] += ww/(abs(iy)+abs(ix)+abs(iz)+1);
//#else
//                u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz] = sf_cadd(u[rtm->sou_y+iy][rtm->sou_x+ix][rtm->sou_z+iz],sf_crmul(ww,1./(abs(iy)+abs(ix)+abs(iz)+1)));
//#endif
//            }
//        }
//    }
}

void bel3d_finalize()
/*< finalize bell taper >*/
{
   if(NULL!=bel3d) { free(**bel3d); free(*bel3d); free(bel3d); bel3d=NULL; }
}


void inject3d(sf_complex ***u,
              sf_complex ***d,
              int tt,
              rtm3d rtm)
/*< inject data into 3d cube >*/
{
    int iy,ix;

#ifdef _OPENMP
#pragma omp parallel for			\
        private(iy,ix)                          \
        shared(rtm,u,d,tt)
#endif
    for    (iy=0; iy<rtm->rec_ny; iy++) {
        for(ix=0; ix<rtm->rec_nx; ix++) {
#ifdef SF_HAS_COMPLEX_H
            u[iy*rtm->rec_jy+rtm->rec_oy][ix*rtm->rec_jx+rtm->rec_ox][rtm->rec_dep] += d[iy][ix][tt];
#else
            u[iy*rtm->rec_jy+rtm->rec_oy][ix*rtm->rec_jx+rtm->rec_ox][rtm->rec_dep] = sf_cadd(u[iy*rtm->rec_jy+rtm->rec_oy][ix*rtm->rec_jx+rtm->rec_ox][rtm->rec_dep],d[iy][ix][tt]);
#endif
        }
    }
}

void extract3d(sf_complex ***u,
               sf_complex ***d,
               int tt,
               rtm3d rtm)
/*< extract data from 3d cube >*/
{
    int iy,ix;

#ifdef _OPENMP
#pragma omp parallel for			\
        private(iy,ix)                          \
        shared(rtm,u,d,tt)
#endif
    for    (iy=0; iy<rtm->rec_ny; iy++) {
        for(ix=0; ix<rtm->rec_nx; ix++) {
            d[iy][ix][tt] = u[iy*rtm->rec_jy+rtm->rec_oy][ix*rtm->rec_jx+rtm->rec_ox][rtm->rec_dep];
        }
    }
}

void mute3d(sf_complex ***d,
            fdm3d fdm,
            rtm3d rtm)
/*< mute direct arrivals for 3d shot record >*/
{
    int iy,ix,it;
    float dist,dist_z,dist_x,dist_y,wt,t1,tt;

    dist_z = fdm->dz*(rtm->rec_dep - rtm->sou_z);
#ifdef _OPENMP
#pragma omp parallel for                              \
        private(iy,ix,it,dist_y,dist_x,dist,t1,tt,wt) \
        shared(rtm,d,dist_z)
#endif
    for(iy=0; iy<rtm->rec_ny; iy++) {
        dist_y = fdm->dy*(iy*rtm->rec_jy+rtm->rec_oy - rtm->sou_y);
        for(ix=0; ix<rtm->rec_nx; ix++) {
            dist_x = fdm->dx*(ix*rtm->rec_jx+rtm->rec_ox - rtm->sou_x);
            dist = sqrtf( powf(dist_z,2.) + powf(dist_x,2.) + powf(dist_y,2.) );
            t1 = rtm->sou_t0 + dist/rtm->vel_w; // first arrival
            for(it=0; it<rtm->rec_nt; it++) {
                tt = it*rtm->rec_jt*rtm->sou_dt;
                if (tt < t1) {
                    wt = expf(1000*(tt-t1));
#ifdef SF_HAS_COMPLEX_H
                    d[iy][ix][it] *= wt;
#else
                    d[iy][ix][it] = sf_crmul(d[iy][ix][it],wt);
#endif
                }
            }
        }
    }
}

void forward(sf_complex ***u,
             fdm3d fdm,
             dft3d dft,
             lrk3d lrk,
             sponge spo)
/*< forward propagate >*/
{
    /* apply lowrank matrices */
    lrk3d_apply(u[0][0], u[0][0], false, fdm, dft, lrk);

    /* apply abc */
    if (NULL!=spo) sponge3d_apply_complex(u, spo, fdm);
}

void reverse(sf_complex ***bu,
             fdm3d fdm,
             dft3d dft,
             lrk3d lrk,
             sponge spo)
/*< reverse and apply imaging >*/
{
    lrk3d_apply(bu[0][0], bu[0][0], true, fdm, dft, lrk);

    /* apply abc */
    if (NULL!=spo) sponge3d_apply_complex(bu, spo, fdm);
}

void ccr(sf_complex ***img,
         sf_complex ***u,
         sf_complex ***bu,
         fdm3d fdm)
/*< cross-correlation imaging condition >*/
{
    int iy,ix,iz;

    if(fdm->ny>1) {
#ifdef _OPENMP
#pragma omp parallel for    \
        private(iy,ix,iz)   \
        shared(img,u,bu,fdm)
#endif
        for (iy=0; iy<fdm->ny; iy++)         {
            for (ix=0; ix<fdm->nx; ix++)     {
                for (iz=0; iz<fdm->nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
                    img[iy][ix][iz] += conjf(u[iy+fdm->nb][ix+fdm->nb][iz+fdm->nb])*bu[iy+fdm->nb][ix+fdm->nb][iz+fdm->nb];
#else
                    img[iy][ix][iz] = sf_cadd(img[iy][ix][iz],sf_cmul(conjf(u[iy+fdm->nb][ix+fdm->nb][iz+fdm->nb]),bu[iy+fdm->nb][ix+fdm->nb][iz+fdm->nb]));
#endif
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for    \
        private(ix,iz)      \
        shared(img,u,bu,fdm)
#endif
        for (ix=0; ix<fdm->nx; ix++)     {
            for (iz=0; iz<fdm->nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
                img[0][ix][iz] += conjf(u[0][ix+fdm->nb][iz+fdm->nb])*bu[0][ix+fdm->nb][iz+fdm->nb];
#else
                img[0][ix][iz] = sf_cadd(img[0][ix][iz],sf_cmul(conjf(u[0][ix+fdm->nb][iz+fdm->nb]),bu[0][ix+fdm->nb][iz+fdm->nb]));
#endif
            }
        }
    }
}

void rtm_finalize()
/*< finalize rtm >*/
{
}

void itoa(int n, char *s)
/*< convert integer to char >*/
{
    int i,j,sign;
    char c;
    sign = n;
    if (n < 0) n=-n;
    i = 0;
    do {
	s[i++] = n%10+'0';
	n = (int) n/10;
    } while(n > 0);
    if (sign <0) s[i++] = '-';
    s[i] = '\0';
    for (i=0,j=strlen(s)-1;i<j;i++,j--){
        c = s[i];
        s[i]=s[j];
        s[j]=c;
    }
}

void setval_complex(sf_complex *u, int n, sf_complex val)
/*< set value >*/
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for                \
    private(i)                          \
    shared(u)
#endif
    for (i=0; i<n; i++) u[i] = val;
}
