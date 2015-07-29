#include "acdpml_gfdm.h"

/*--- time step functions ---------------------------------------------------*/
extern "C"{
    void acdpml_2d_2_d(float **uc,   float **ucd,
                       float **up,   float **upd,
                       float **csq,  float **csqd,
                       float **phi1, float **phi1d,
                       float **phi0, float **phi0d,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0, float *c1,
                       int *lbc, int *rbc);
    void acdpml_2d_4_d(float **uc,   float **ucd,
                       float **up,   float **upd,
                       float **csq,  float **csqd,
                       float **phi1, float **phi1d,
                       float **phi0, float **phi0d,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0, float *c1, float *c2,
                       int *lbc, int *rbc);
    void acdpml_2d_8_d(float **uc,   float **ucd,
                       float **up,   float **upd,
                       float **csq,  float **csqd,
                       float **phi1, float **phi1d,
                       float **phi0, float **phi0d,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0,  float *c1, float *c2,
                       float *c3, float *c4,
                       int *lbc, int *rbc);
    
    void acdpml_2d_2_b(float **uc,   float **ucb,
                       float **up,   float **upb,
                       float **csq,  float **csqb,
                       float **phi1, float **phi1b,
                       float **phi0, float **phi0b,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0, float *c1,
                       int *lbc, int *rbc);
    void acdpml_2d_4_b(float **uc,   float **ucb,
                       float **up,   float **upb,
                       float **csq,  float **csqb,
                       float **phi1, float **phi1b,
                       float **phi0, float **phi0b,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0, float *c1, float *c2,
                       int *lbc, int *rbc);
    void acdpml_2d_8_b(float **uc,   float **ucb,
                       float **up,   float **upb,
                       float **csq,  float **csqb,
                       float **phi1, float **phi1b,
                       float **phi0, float **phi0b,
                       float *dp1,   float *dp0,
                       float *di,    float dt,
                       int *s,       int *e,
                       float c0,  float *c1, float *c2,
                       float *c3, float *c4,
                       int *lbc, int *rbc);
}


int acdpml_tsfm(RDOM * p, RDOM * r, int ia, void * fdpars) {
    
    // swap pointer workspace
    ireal tmp;
    IPNT i;
    
    // pointers for 2D case
    register ireal ** restrict uc2;
    register ireal ** restrict up2;
    register ireal ** restrict csq2;
    register ireal ** restrict phi1;
    register ireal ** restrict phi0;
    
    register ireal ** restrict uc2d;
    register ireal ** restrict up2d;
    register ireal ** restrict csq2d;
    register ireal ** restrict phi1d;
    register ireal ** restrict phi0d;
    
    /* pointers for 3D case
     register ireal *** restrict uc3;
     register ireal *** restrict up3;
     register ireal *** restrict csq3;
     register ireal *** restrict uc3d;
     register ireal *** restrict up3d;
     register ireal *** restrict csq3d;
     */
    
    int ndim;                       // problem dmn
    IPNT s, s0;                     // loop starts
    IPNT e, e0;                     // loop ends
    
    // acd struct
    ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)fdpars;
    
    // extract dimn info
    ra_ndim(&(r->_s[D_UC]),&ndim);
    ra_gse(&(r->_s[D_UC]),s,e);
    ra_a_gse(&(r->_s[D_UC]),s0,e0);
    
    if (ndim == 2) {
        
        // 2D computational arrays
        uc2    = (r->_s)[D_UC ]._s2;
        up2    = (r->_s)[D_UP ]._s2;
        csq2   = (r->_s)[D_CSQ]._s2;
        phi1  = (r->_s)[D_PHI1]._s2;
        phi0  = (r->_s)[D_PHI0]._s2;
        
        uc2d   = (p->_s)[D_UC ]._s2;
        up2d   = (p->_s)[D_UP ]._s2;
        csq2d  = (p->_s)[D_CSQ]._s2;
        phi1d  = (p->_s)[D_PHI1]._s2;
        phi0d  = (p->_s)[D_PHI0]._s2;
        
        // 2nd order case
        if (acdpmlpars->k == 1) {
            acdpml_2d_2_d(uc2,  uc2d,
             up2,  up2d,
             csq2, csq2d,
             phi1, phi1d,
             phi0, phi0d,
             acdpmlpars->dp1,
             acdpmlpars->dp0,
             acdpmlpars->di,
             acdpmlpars->dt,
             s, e,
             acdpmlpars->c0,
             acdpmlpars->c1,
             acdpmlpars->lbc,
             acdpmlpars->rbc);
                 }
        else if (acdpmlpars->k == 2) {
            acdpml_2d_4_d(uc2,  uc2d,
                          up2,  up2d,
                          csq2, csq2d,
                          phi1, phi1d,
                          phi0, phi0d,
                          acdpmlpars->dp1,
                          acdpmlpars->dp0,
                          acdpmlpars->di,
                          acdpmlpars->dt,
                          s, e,
                          acdpmlpars->c0,
                          acdpmlpars->c1,
                          acdpmlpars->c2,
                          acdpmlpars->lbc,
                          acdpmlpars->rbc);
        }
        else if (acdpmlpars->k == 4) {
            acdpml_2d_8_d(uc2,  uc2d,
                          up2,  up2d,
                          csq2, csq2d,
                          phi1, phi1d,
                          phi0, phi0d,
                          acdpmlpars->dp1,
                          acdpmlpars->dp0,
                          acdpmlpars->di,
                          acdpmlpars->dt,
                          s, e,
                          acdpmlpars->c0,
                          acdpmlpars->c1,
                          acdpmlpars->c2,
                          acdpmlpars->c3,
                          acdpmlpars->c4,
                          acdpmlpars->lbc,
                          acdpmlpars->rbc);
        }
        else {
            fprintf(stderr,"ERROR: acdpml_tsfm\n");
            fprintf(stderr,"called with half-order != 1\n");
            return E_BADINPUT;
        }
        
        for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
            for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
                tmp=uc2d[i[1]][i[0]];
                uc2d[i[1]][i[0]]=up2d[i[1]][i[0]];
                up2d[i[1]][i[0]]=tmp;
                tmp=uc2[i[1]][i[0]];
                uc2[i[1]][i[0]]=up2[i[1]][i[0]];
                up2[i[1]][i[0]]=tmp;
            }
        }
    }
    else if (ndim == 3) {
        
        /*
         uc3    = (r->_s)[D_UC ]._s3;
         up3    = (r->_s)[D_UP ]._s3;
         csq3   = (r->_s)[D_CSQ]._s3;
         uc3d   = (p->_s)[D_UC ]._s3;
         up3d   = (p->_s)[D_UP ]._s3;
         csq3d  = (p->_s)[D_CSQ]._s3;
         */
        
        // 2nd order case
        if (acdpmlpars->k == 1) {
            /*
             acd_3d_2(uc3, up3, csq3,
             s, e,
             acdpars->c0,
             acdpars->c1);
             */
        }
        // 4th order case
        else if (acdpmlpars->k == 2) {
            /*
             acd_3d_4(uc3, up3, csq3,
             s, e,
             acdpars->c0,
             acdpars->c1, acdpars->c2,
             acdpars->lbc, acdpars->rbc);
             */
        }
        // 8th order case
        else if (acdpmlpars->k == 4) {
            /*
             acd_3d_8(uc3, up3, csq3,
             s, e,
             acdpars->c0,
             acdpars->c1, acdpars->c2,
             acdpars->c3, acdpars->c4,
             acdpars->lbc, acdpars->rbc);
             */
        }
        else {
            fprintf(stderr,"ERROR: acdpml_tsfm\n");
            fprintf(stderr,"called with half-order != 1\n");
            return E_BADINPUT;
        }
        
        for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
            for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
                for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
                    tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
                    ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
                    ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
                }
            }
        }
    }
    else {
        fprintf(stderr,"ERROR: acd_step\n");
        fprintf(stderr,"called with space dim != 2 or 3\n");
        return E_BADINPUT;
    }
    
    return 0;
}

int acdpml_tsam(RDOM * p, RDOM * r, int ia, void * fdpars) {
    
    // swap pointer workspace
    ireal tmp;
    IPNT i;
    IPNT n;
    
    // pointers for 2D case
    register ireal ** restrict uc2;
    register ireal ** restrict up2;
    register ireal ** restrict csq2;
    register ireal ** restrict phi1;
    register ireal ** restrict phi0;
    
    register ireal ** restrict uc2b;
    register ireal ** restrict up2b;
    register ireal ** restrict csq2b;
    register ireal ** restrict phi1b;
    register ireal ** restrict phi0b;
    
    /* pointers for 3D case
     register ireal *** restrict uc3;
     register ireal *** restrict up3;
     register ireal *** restrict csq3;
     register ireal *** restrict uc3d;
     register ireal *** restrict up3d;
     register ireal *** restrict csq3d;
     */
    
    int ndim;                       // problem dmn
    IPNT s, s0;                     // loop starts
    IPNT e, e0;                     // loop ends
    // acd struct
    ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)fdpars;
    
    // extract dimn info
    ra_ndim(&(r->_s[D_UC]),&ndim);
    ra_gse(&(r->_s[D_UC]),s,e);
    ra_a_gse(&(r->_s[D_UC]),s0,e0);
    ra_a_size(&(r->_s[D_UC]),n);
    
    if (ndim == 2) {
        // 2D computational arrays
        uc2    = (r->_s)[D_UC ]._s2;
        up2    = (r->_s)[D_UP ]._s2;
        csq2   = (r->_s)[D_CSQ]._s2;
        phi1  = (r->_s)[D_PHI1]._s2;
        phi0  = (r->_s)[D_PHI0]._s2;
        uc2b   = (p->_s)[D_UC ]._s2;
        up2b   = (p->_s)[D_UP ]._s2;
        csq2b  = (p->_s)[D_CSQ]._s2;
        phi1b  = (p->_s)[D_PHI1]._s2;
        phi0b  = (p->_s)[D_PHI0]._s2;
        for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
            for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
                tmp=uc2b[i[1]][i[0]];
                uc2b[i[1]][i[0]]=up2b[i[1]][i[0]];
                up2b[i[1]][i[0]]=tmp;
            }
        }
        
        // 2nd order case
        if (acdpmlpars->k == 1) {
            acdpml_2d_2_b(uc2,  uc2b,
             up2,  up2b,
             csq2, csq2b,
             phi1, phi1b,
             phi0, phi0b,
             acdpmlpars->dp1,
             acdpmlpars->dp0,
             acdpmlpars->di,
             acdpmlpars->dt,
             s, e,
             acdpmlpars->c0,
             acdpmlpars->c1,
             acdpmlpars->lbc,
             acdpmlpars->rbc);
        }
        else if (acdpmlpars->k == 2){
            acdpml_2d_4_b(uc2,  uc2b,
                          up2,  up2b,
                          csq2, csq2b,
                          phi1, phi1b,
                          phi0, phi0b,
                          acdpmlpars->dp1,
                          acdpmlpars->dp0,
                          acdpmlpars->di,
                          acdpmlpars->dt,
                          s, e,
                          acdpmlpars->c0,
                          acdpmlpars->c1,
                          acdpmlpars->c2,
                          acdpmlpars->lbc,
                          acdpmlpars->rbc);
        }
        else if (acdpmlpars->k == 4){
            acdpml_2d_8_b(uc2,  uc2b,
                          up2,  up2b,
                          csq2, csq2b,
                          phi1, phi1b,
                          phi0, phi0b,
                          acdpmlpars->dp1,
                          acdpmlpars->dp0,
                          acdpmlpars->di,
                          acdpmlpars->dt,
                          s, e,
                          acdpmlpars->c0,
                          acdpmlpars->c1,
                          acdpmlpars->c2,
                          acdpmlpars->c3,
                          acdpmlpars->c4,
                          acdpmlpars->lbc,
                          acdpmlpars->rbc);
        }
        else {
            fprintf(stderr,"ERROR: acdpml_tsam\n");
            fprintf(stderr,"called with half-order != 1\n");
            return E_BADINPUT;
        }
    }
    
    else if (ndim == 3) {
        
        for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
            for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
                for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
                    tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
                    ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
                    ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
                }
            }
        }
        
        /*
         uc3    = (r->_s)[D_UC ]._s3;
         up3    = (r->_s)[D_UP ]._s3;
         csq3   = (r->_s)[D_CSQ]._s3;
         uc3d   = (p->_s)[D_UC ]._s3;
         up3d   = (p->_s)[D_UP ]._s3;
         csq3d  = (p->_s)[D_CSQ]._s3;
         */
        
        // 2nd order case 
        if (acdpmlpars->k == 1) {
            /*
             acd_3d_2(uc3, up3, csq3, 
             s, e, 
             acdpars->c0, 
             acdpars->c1);
             */
        }
        // 4th order case
        else if (acdpmlpars->k == 2) {
            /*
             acd_3d_4(uc3, up3, csq3,
             s, e, 
             acdpars->c0, 
             acdpars->c1, acdpars->c2,
             acdpars->lbc, acdpars->rbc);
             */
        }
        // 8th order case
        else if (acdpmlpars->k == 4) {
            /*
             acd_3d_8(uc3, up3, csq3,
             s, e, 
             acdpars->c0, 
             acdpars->c1, acdpars->c2,
             acdpars->c3, acdpars->c4,
             acdpars->lbc, acdpars->rbc);
             */
        } 
        else {
            fprintf(stderr,"ERROR: acdpml_tsam\n");
            fprintf(stderr,"called with half-order != 1\n");
            return E_BADINPUT;
        }
        
    }
    else {
        fprintf(stderr,"ERROR: acdpml_tsam\n");
        fprintf(stderr,"called with space dim != 2 or 3\n");
        return E_BADINPUT;
    }
    
    return 0;
}


