#include "acd_gfdm.h"

/*--- time step functions ---------------------------------------------------*/
extern "C" {
void acd_2d_2_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1);

void acd_2d_2_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1);

void acd_2d_4_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1, float *c2, 
		int *lbc, int *rbc);

void acd_2d_4_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1, float *c2, 
		int *lbc, int *rbc);

void acd_2d_8_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1, float *c2, 
		float *c3, float *c4, int *lbc, int *rbc);

void acd_2d_8_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1, float *c2, 
		float *c3, float *c4, int *lbc, int *rbc);
}


int acd_tsfm(RDOM * p, RDOM * r, int ia, void * fdpars) {

  // swap pointer workspace
  ireal tmp;
  IPNT i;

  // pointers for 2D case
  register ireal ** restrict uc2;
  register ireal ** restrict up2;
  register ireal ** restrict csq2;
  register ireal ** restrict uc2d;
  register ireal ** restrict up2d;
  register ireal ** restrict csq2d;

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
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;

  // extract dimn info
  ra_ndim(&(r->_s[D_UC]),&ndim);
  ra_gse(&(r->_s[D_UC]),s,e);
  ra_a_gse(&(r->_s[D_UC]),s0,e0);

  if (ndim == 2) {

    // 2D computational arrays
    uc2    = (r->_s)[D_UC ]._s2;
    up2    = (r->_s)[D_UP ]._s2;
    csq2   = (r->_s)[D_CSQ]._s2;
    uc2d   = (p->_s)[D_UC ]._s2;
    up2d   = (p->_s)[D_UP ]._s2;
    csq2d  = (p->_s)[D_CSQ]._s2;

    // 2nd order case 
    if (acdpars->k == 1) {
      acd_2d_2_d(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d, 
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1);
    }
    // 4th order case
    else if (acdpars->k == 2) {
      acd_2d_4_d(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->lbc, acdpars->rbc);
    }
    // 8th order case
    else if (acdpars->k == 4) {
      acd_2d_8_d(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->c3, acdpars->c4,
		 acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
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
    if (acdpars->k == 1) {
      /*
      acd_3d_2(uc3, up3, csq3, 
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1);
      */
    }
    // 4th order case
    else if (acdpars->k == 2) {
      /*
      acd_3d_4(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->lbc, acdpars->rbc);
      */
    }
    // 8th order case
    else if (acdpars->k == 4) {
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
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
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

int acd_tsam(RDOM * p, RDOM * r, int ia, void * fdpars) {

  // swap pointer workspace
  ireal tmp;
  IPNT i;
  IPNT n;

  // pointers for 2D case
  register ireal ** restrict uc2;
  register ireal ** restrict up2;
  register ireal ** restrict csq2;
  register ireal ** restrict uc2d;
  register ireal ** restrict up2d;
  register ireal ** restrict csq2d;

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
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;

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
    uc2d   = (p->_s)[D_UC ]._s2;
    up2d   = (p->_s)[D_UP ]._s2;
    csq2d  = (p->_s)[D_CSQ]._s2;
    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	tmp=uc2d[i[1]][i[0]];
	uc2d[i[1]][i[0]]=up2d[i[1]][i[0]];
	up2d[i[1]][i[0]]=tmp;
      }
    }

    // 2nd order case 
    if (acdpars->k == 1) {
      acd_2d_2_b(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d, 
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1);
    }
    // 4th order case
    else if (acdpars->k == 2) {
      acd_2d_4_b(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->lbc, acdpars->rbc);
    }
    // 8th order case
    else if (acdpars->k == 4) {
      acd_2d_8_b(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->c3, acdpars->c4,
		 acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
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
    if (acdpars->k == 1) {
      /*
      acd_3d_2(uc3, up3, csq3, 
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1);
      */
    }
    // 4th order case
    else if (acdpars->k == 2) {
      /*
      acd_3d_4(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->lbc, acdpars->rbc);
      */
    }
    // 8th order case
    else if (acdpars->k == 4) {
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
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

  }
  else {
    fprintf(stderr,"ERROR: acd_step\n");
    fprintf(stderr,"called with space dim != 2 or 3\n");
    return E_BADINPUT;
  }
  
  return 0;
}


