/* Time reversal imaging of passive seismic data linear operator */
#include <rsf.h>
#include "triutil2.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#ifndef _triutil2_h

#define NOP 4 /* derivative operator half-size */
#define C0 -205.0f/72.0f
#define C1 +8.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 +8.0f/315.0f
#define C4 -1.0f/560.0f
#define Lap(a,ix,iz,sx,sz,v)  ( ( C4*(a[ix+4][iz  ] + a[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ] + a[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ] + a[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ] + a[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]) )*sx            +      \
                                ( C4*(a[ix  ][iz+4] + a[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3] + a[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2] + a[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1] + a[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]) )*sz )*v[ix][iz]
#define LapT(a,ix,iz,sx,sz,v) ( ( C4*(a[ix+4][iz  ]*v[ix+4][iz  ] + a[ix-4][iz  ]*v[ix-4][iz  ]) +      \
                                  C3*(a[ix+3][iz  ]*v[ix+3][iz  ] + a[ix-3][iz  ]*v[ix-3][iz  ]) +      \
                                  C2*(a[ix+2][iz  ]*v[ix+2][iz  ] + a[ix-2][iz  ]*v[ix-2][iz  ]) +      \
                                  C1*(a[ix+1][iz  ]*v[ix+1][iz  ] + a[ix-1][iz  ]*v[ix-1][iz  ]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sx                          +      \
                                ( C4*(a[ix  ][iz+4]*v[ix  ][iz+4] + a[ix  ][iz-4]*v[ix  ][iz-4]) +      \
                                  C3*(a[ix  ][iz+3]*v[ix  ][iz+3] + a[ix  ][iz-3]*v[ix  ][iz-3]) +      \
                                  C2*(a[ix  ][iz+2]*v[ix  ][iz+2] + a[ix  ][iz-2]*v[ix  ][iz-2]) +      \
                                  C1*(a[ix  ][iz+1]*v[ix  ][iz+1] + a[ix  ][iz-1]*v[ix  ][iz-1]) +      \
                                  C0*(a[ix  ][iz  ]*v[ix  ][iz  ]) )*sz )
/*^*/

typedef struct tri2 *tri2d;
/*^*/

struct tri2{
    bool verb,abc;
    int  nt, nx, nz, nb, depth, nxpad, nzpad, nzxpad, ngeo, ngrp, **geo;
    float dt2, idz2, idx2, cb;
};
/*^*/

#endif

/***************************************************************/
tri2d tri2d_make(bool verb, bool abc,
                int nt, int nx, int nz, int nb, int depth, int ngeo, int ngrp, int **geo,
                float dt, float dx, float dz, float cb)
/*< initialize tri2d utilities >*/
{   
    tri2d tri;
    int ig;

    tri = (tri2d) sf_alloc(1,sizeof(*tri));

    tri->verb  = verb;
    tri->abc   = abc;
    tri->nt    = nt;
    tri->nx    = nx;
    tri->nz    = nz;
    tri->nb    = nb;
    tri->depth = depth;
    tri->dt2   = dt*dt;
    tri->idx2  = 1.0f/(dx*dx);
    tri->idz2  = 1.0f/(dz*dz);
    tri->cb    = cb; 
    tri->ngeo  = ngeo;
    tri->ngrp  = ngrp;
    tri->geo   = geo;

    tri->nxpad = nx+2*nb;
    tri->nzpad = nz+2*nb;
    tri->nzxpad= tri->nzpad*tri->nxpad;
    tri->depth = depth+nb;
    if (NULL!=tri->geo) {
        for (ig=0; ig<tri->ngeo; ig++) {
            tri->geo[ig][2] += nb;
            tri->geo[ig][1] += nb;
        }
    }

    return tri;
}

/***************************************************************/
static tri2d tr;
static float **vvpad, **u0, **u1, **u2, **tmp;

void timerev_init(bool verb, bool abc,
                  int nt, int nx, int nz, int nb, int depth, int ngeo, int ngrp, int **geo,
                  float dt, float dx, float dz, float cb,
                  float **vv)
/*< initialize >*/
{
    int ix, iz;

    tr = tri2d_make(verb, abc, nt, nx, nz, nb, depth, ngeo, ngrp, geo, dt, dx, dz, cb);

    /* set Laplacian coefficients */
    vvpad = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u0    = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u1    = sf_floatalloc2(tr->nzpad, tr->nxpad);
    u2    = sf_floatalloc2(tr->nzpad, tr->nxpad);

    /* pad boundary */
    for     (ix=0; ix<tr->nx; ix++)
        for (iz=0; iz<tr->nz; iz++)
            vvpad[ix+tr->nb][iz+tr->nb] = vv[ix][iz]*vv[ix][iz]*tr->dt2;
    for     (ix=0; ix<tr->nxpad; ix++){
        for (iz=0; iz<tr->nb;    iz++){
            vvpad[ix][          iz  ] = vvpad[ix][          tr->nb  ];
            vvpad[ix][tr->nzpad-iz-1] = vvpad[ix][tr->nzpad-tr->nb-1];
        }
    }
    for     (ix=0; ix<tr->nb;    ix++){
        for (iz=0; iz<tr->nzpad; iz++){
            vvpad[          ix  ][iz]=vvpad[          tr->nb  ][iz];
            vvpad[tr->nxpad-ix-1][iz]=vvpad[tr->nxpad-tr->nb-1][iz];
        }
    }

    /* absorbing boundary condition */
    abc_init(tr);
}

void inject(int igrp, int it, float **uu, float **dd)
/*< inject data >*/
{
    int igeo;
    /* loop over geophones */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(igeo)
#endif
    for (igeo=0; igeo<tr->ngeo; igeo++) {
        if ( tr->geo[igeo][0] == igrp )
            uu[tr->geo[igeo][2]][tr->geo[igeo][1]] += dd[igeo][it];
    }
}

void extract(int igrp, int it, float **uu, float **dd)
/*< extract data >*/
{
    int igeo;
    /* loop over geophones */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(igeo)
#endif
    for (igeo=0; igeo<tr->ngeo; igeo++) {
        if ( tr->geo[igeo][0] == igrp )
            dd[igeo][it] = uu[tr->geo[igeo][2]][tr->geo[igeo][1]];
    }
}

void timerev_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< time reversal imaging linear operator >*/
{
    int ix, iz, it, ig;
    float **dd, ***ww;

    if (nm!=tr->nz*tr->nx*tr->nt || nd!=tr->nt*tr->ngeo) sf_error("%s: wrong dimensions",__FILE__);
    sf_adjnull(adj, add, nm, nd, mod, dat);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nxpad; ix++)
                for (iz=0; iz<tr->nzpad; iz++)
                {
                    u0[ix][iz] = 0.0f;
                    u1[ix][iz] = 0.0f;
                    u2[ix][iz] = 0.0f;
                }
 
    /* map 1d to 2d */
    dd = (float**) sf_alloc (tr->ngeo,sizeof(float*)); 
    dd[0] = dat;
    for (ix=1; ix<tr->ngeo; ix++) dd[ix] = dd[0]+ix*tr->nt;

    /* map 1d to 3d */
    ww = (float***) sf_alloc (tr->nt,sizeof(float**));
    ww[0] = (float**) sf_alloc (tr->nx*tr->nt,sizeof(float*));
    ww[0][0] = mod;
    for (ix=1; ix<tr->nx*tr->nt; ix++) ww[0][ix] = ww[0][0]+ix*tr->nz; 
    for (it=1; it<tr->nt; it++) ww[it] = ww[0]+it*tr->nx;

    if (adj) { /* migration */
        
        for (it=tr->nt-1; it>-1; it--){
           
            /* 4 - apply abc */
            if (tr->abc) abc_apply(u1[0],tr);
            if (tr->abc) abc_apply(u0[0],tr);

            /* 3 - image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    ww[it][ix][iz] += u1[ix+tr->nb][iz+tr->nb]*vvpad[ix+tr->nb][iz+tr->nb];

            /* 2 - time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = LapT(u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;

            /* 1 - inject data */
            if (NULL!=tr->geo) {
                for (ig=0; ig<tr->ngrp; ig++)
                    inject(ig, it, u1, dd);
            } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
                for (ix=tr->nb; ix<tr->nb+tr->nx; ix++)
                    u1[ix][tr->depth] += dd[ix-tr->nb][it];
            }
 
        } /* it loop */

    } else { /* modeling */
    	
    	for (it=0; it<tr->nt; it++){

             /* 1 - record data */
            if (NULL!=tr->geo) {
                for (ig=0; ig<tr->ngrp; ig++)
                    extract(ig, it, u1, dd);
            } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
                for (ix=tr->nb; ix<tr->nb+tr->nx; ix++)
                    dd[ix-tr->nb][it] += u1[ix][tr->depth];
            }
           
            /* 2 - time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = Lap(u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;

            /* 3 - inject source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    u1[ix+tr->nb][iz+tr->nb] += ww[it][ix][iz]*vvpad[ix+tr->nb][iz+tr->nb];

            /* 4 - apply abc */
            if (tr->abc) abc_apply(u0[0],tr);
            if (tr->abc) abc_apply(u1[0],tr);

        } /* it loop */
        
    }

    free(dd);
    free(*ww); free(ww);
}

void timerev_close()
/*< finalize >*/
{
    abc_close();
    free(tr);
    free(*vvpad); free(vvpad);
    free(*u0); free(u0);
    free(*u1); free(u1);
    free(*u2); free(u2);
}

/***************************************************************/
void ctimerev(int ngrp, float ***ww, float **dd)
/*< correlative time reversal imaging condition >*/
{
    int ix, iz, it, ig, counter, *beg=NULL, *end=NULL;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
    for         (it=0; it<tr->nt; it++)
        for     (ix=0; ix<tr->nx; ix++)
            for (iz=0; iz<tr->nz; iz++)
                ww[it][ix][iz] = 1.0f;

    if (NULL==tr->geo) {
        /* set start and end index */
        beg=sf_intalloc(ngrp);
        end=sf_intalloc(ngrp);
        counter = 0;
        for (ig=0; ig<ngrp; ig++) {
            beg[ig] = counter;
            counter += tr->nx/ngrp;
            end[ig] = counter;
        }
        end[ngrp-1] = tr->nx;
        if (tr->verb) {
            for (ig=0; ig<ngrp; ig++) {
                sf_warning("beg[%d]=%d",ig,beg[ig]);
                sf_warning("end[%d]=%d",ig,end[ig]);
            }
        }
    }

    for (ig=0; ig<ngrp; ig++) { /* loop over subgroups of receivers */

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nxpad; ix++)
                for (iz=0; iz<tr->nzpad; iz++)
                {
                    u0[ix][iz] = 0.0f;
                    u1[ix][iz] = 0.0f;
                    u2[ix][iz] = 0.0f;
                }

        for (it=tr->nt-1; it>-1; it--){
            if (tr->verb) sf_warning("Time reversal: %d/%d;", it, 0);

            /* time stepping */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=NOP; ix<tr->nxpad-NOP; ix++){
                for (iz=NOP; iz<tr->nzpad-NOP; iz++){
                    u2[ix][iz] = Lap (u1,ix,iz,tr->idx2,tr->idz2,vvpad) + 2.0f*u1[ix][iz] - u0[ix][iz];
                }
            }
            /* rotate pointers */
            tmp=u0; u0=u1; u1=u2; u2=tmp;
            if (tr->abc) abc_apply(u1[0],tr);
            if (tr->abc) abc_apply(u0[0],tr);

            /* inject data */
            if (NULL!=tr->geo) {
                inject(ig, it, u1, dd);
            } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
                for (ix=tr->nb+beg[ig]; ix<tr->nb+end[ig]; ix++)
                    u1[ix][tr->depth] += dd[ix-tr->nb][it];
            }

            /* image source */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix, iz)
#endif
            for     (ix=0; ix<tr->nx; ix++)
                for (iz=0; iz<tr->nz; iz++)
                    ww[it][ix][iz] *= u1[ix+tr->nb][iz+tr->nb];

        } /* it loop */
        if (tr->verb) sf_warning(".");

    } /* ig loop */


}

/***************************************************************/
/* absorbing boundary */
static float *decay=NULL;

void abc_init(tri2d tri)
/*< initialization >*/
{
    if(tri->nb) decay =  sf_floatalloc(tri->nb);
    abc_cal(tri->nb,tri->cb,decay);
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(NULL!=decay) free(decay);
    decay = NULL;
}

void abc_apply(float *a /*2-D matrix*/,
               tri2d tri) 
/*< boundary decay>*/
{
    int iz, ix;

#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix)
{
#endif
    /* top & bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < tri->nb; iz++) {  
        for (ix=0; ix < tri->nxpad; ix++) {
	  a[tri->nzpad*ix +              iz] *= decay[iz];
	  a[tri->nzpad*ix + tri->nzpad-1-iz] *= decay[iz];
        }
    }
    /* left & right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < tri->nzpad; iz++) {  
        for (ix=0; ix < tri->nb; ix++) {
	  a[tri->nzpad*              ix  + iz] *= decay[ix];
	  a[tri->nzpad*(tri->nxpad-1-ix) + iz] *= decay[ix];
        }
    }
#ifdef _OPENMP
}
#endif
}

void abc_cal(int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< calculate absorbing coefficients >*/
{
    int ib;
    if(!nb) return;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
    for(ib=0; ib<nb; ib++){
        w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
    }
}

/***************************************************************/
/* sliding window normalization */

void threshold(bool step, int n, float hard, float *dat)
/*< in-place hard thresholding >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++) {
        if (dat[i]<hard) {
            if (step) dat[i] = 0.0f;
            else { 
                /*dat[i] /= hard;*/
                dat[i] = 0.5*(1.+cosf(SF_PI*(dat[i]/hard-1.))); /* Tukey window */
            }
        } else dat[i] = 1.0f;
    }
}

void absval(int n, float *dat)
/*< in-place absolute value >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++)
        dat[i] = fabs(dat[i]);
}

void autopow(int n, float p, float *dat)
/*< in-place auto-correlation with abs >*/
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<n; i++)
        dat[i] = powf(fabs(dat[i]),p);
}

float maxval(int n, float *dat)
/*< maximum absolute value and variance (optional) >*/
{
    /* no parallelism */
    float dd, max=0;
    int i;

    for (i=0; i<n; i++) { 
        dd = dat[i];
        if (max<dd) max = dd;
    }

    return max;
}

void scale(float a, int n, float *dat)
/*< scale an array >*/
{
    /* no parallelism */
    int i;

    for (i=0; i<n; i++)
        dat[i] *= a;
}

void swnorm(bool verb, bool sw, int nz, int nx, int nt, int size, float perc, float *dat)
/*< local (sliding-window) normalization >*/
{
    int i, nzx, nzxt;
    float *dat0,den,factor;
    float max_all,pad;

    nzx = nz*nx;
    nzxt = nzx*nt;

    max_all = maxval(nzxt,dat);
    pad = max_all*perc/100.0f;
    if(verb) sf_warning("max_all=%g",max_all);

    dat0 = sf_floatalloc(nzxt);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<nzxt; i++) dat0[i] = dat[i];

    if (!sw) {

        scale(1./max_all,nzxt,dat);

    } else {

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=0; i<size; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(i+1+size),dat0);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=size; i<nt-size; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(2*size+1),dat0+(i-size)*nzx);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,den,factor)
#endif
        for (i=nt-size; i<nt; i++) {
            if (verb) sf_warning("i = %d/%d;",i,nt);
            den = maxval(nzx*(size+1+(nt-1-i)),dat0+(i-size)*nzx);
            if (den <= pad) den = pad;
            factor = 1.0f/den;
            scale(factor,nzx,dat+i*nzx);
        }
        if (verb) sf_warning(".");

    }

    free(dat0);
}

/*********************************************************/
/* smoothing */
void smooth(int n1, int n2, int n3,
            int rect1, int rect2, int rect3,
            int nrep,
	    float *dat)
/*< Generate reflectivity map with smoothing >*/
{   
    int i, j, i0, irep, n123; 
    sf_triangle trig;
    int n[3],s[3],rect[3];

    n[0]=n1; n[1]=n2; n[2]=n3;
    s[0]=1;  s[1]=n1; s[2]=n1*n2;
    rect[0]=rect1; rect[1]=rect2; rect[2]=rect3;
    n123=n1*n2*n3;
    
    /* 2-d triangle smoothing */
    for (i=0;i<3;i++) {
        if (rect[i] <= 1) continue;
        trig = sf_triangle_init (rect[i],n[i],false);
        for (j=0; j < n123/n[i]; j++) {
            i0 = sf_first_index (i,j,3,n,s);
            for (irep=0; irep < nrep; irep++) {
                sf_smooth2 (trig,i0,s[i],false,dat);
            }
        }
        sf_triangle_close(trig);
    }
}
