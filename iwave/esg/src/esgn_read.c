/* interface to sg_readtimegrid */

#include "esgn_read.h"
#include "fd.h"
#include "esgn.h"

/**
 * load geometry grid from p-wave, or s-wave velocity, or Lame constants (lambda or mu)
 */

int init_elastic_geom_par(grid * g, PARARRAY par, FILE * stream) {
  char * pvelkey   = NULL;
  char * svelkey   = NULL;
  char * lambdakey = NULL;
  char * mukey     = NULL;
  int err = 0;
  
  int i;
  int tmp;
  char key[3];
  size_t kl = 3;
  
  ps_ffcstring(par, "pvelocity", &pvelkey);
  if (pvelkey) {
    err = read_grid(g, pvelkey, stream);
    free(pvelkey);
    return err;
  }
  
  ps_ffcstring(par, "svelocity", &svelkey);
  if (svelkey) {
    err = read_grid(g, svelkey, stream);
    free(svelkey);
    return err;
  }
  
  ps_ffcstring(par, "lambda", &lambdakey);
  if (lambdakey) {
    err = read_grid(g, lambdakey, stream);
    free(lambdakey);
    return err;
  }

  ps_ffcstring(par, "mu", &mukey);
  if (mukey) {
    err = read_grid(g, mukey, stream);
    free(mukey);
    return err;
  }
  
  /*< otherwise look for n's, d's and o's in parameter list */
  g->dim = 0;
  for (i = 0;i < RARR_MAX_NDIM;i ++) {
    snprintf(key, kl, "n%d", i+1);
    tmp = 1;
    ps_ffint(par, key, &tmp);
    g->axes[i].n = tmp;
    snprintf(key, kl, "d%d", i+1);
    if (ps_ffreal(par, key, &(g->axes[i].d)))  g->axes[i].d=1.0;
    snprintf(key, kl, "o%d", i+1);
    if (ps_ffreal(par, key, &(g->axes[i].o)))  g->axes[i].o=0.0;
    /* determine dim by finding least axis index with n>1 */
    if (g->axes[i].n > 1)  g->dim = iwave_max(g->dim,i+1);
  }
  
  if (RARR_MAX_NDIM > 0) {
    tmp = 1;
    ps_ffint(par, "z_axis", &tmp);
    tmp--;
    if (tmp < 0 || tmp > RARR_MAX_NDIM-1) {
      fprintf(stream, "Error: init_acoustic_geom\n");
      fprintf(stream, "z_axis index = %d out of range for dim = %d\n",
              tmp, RARR_MAX_NDIM);
      return E_OTHER;
    }
    g->axes[tmp].id = 0;
  }

  if (RARR_MAX_NDIM > 1) {
    tmp = 2;
    ps_ffint(par, "x_axis", &tmp);
    tmp--;
    if (tmp < 0 || tmp > RARR_MAX_NDIM-1) {
      fprintf(stream, "Error: init_acoustic_geom\n");
      fprintf(stream, "x_axis index = %d out of range for dim = %d\n",
              tmp, RARR_MAX_NDIM);
      return E_OTHER;
    }
    g->axes[tmp].id = 1;
  }
  
  if (RARR_MAX_NDIM > 2) {
    tmp = 3;
    ps_ffint(par, "y_axis", &tmp);
    tmp--;
    if (tmp < 0 || tmp > RARR_MAX_NDIM-1) {
      fprintf(stream, "Error: init_acoustic_geom\n");
      fprintf(stream, "y_axis index = %d out of range for dim = %d\n",
              tmp, RARR_MAX_NDIM);
      return E_OTHER;
    }
    g->axes[tmp].id = 2;
  }
  return err;
}

/** 
 * compute dt only based on cfl value and maximum velocity
 * provided in para file 
*/
int esg_readtimegrid(PARARRAY * par, FILE * stream, IMODEL * model) {
  
  int i;                        /* counter */
  int ndim;                     /* model dimn */
  int err = 0;                  /* error flag */
  ireal cmax = CMAX_DEF;        /* from input */
  ireal cfl = CFL_DEF;          /* from input */
  ireal a;              
  char * mode = NULL;
  grid g = model->g;
  ndim = g.dim;
  
  ps_ffcstring(*par,"mode",&mode);      
  a = g.axes[0].d;
  
  if (mode==NULL)
    for ( i = 1; i < ndim; ++i ) a = iwave_min(a, g.axes[i].d);
  
  if (ps_ffreal(*par,"cmax",&cmax)) 
    fprintf(stream,"NOTE: esg_readtimegrid - using default cmax = %e\n",cmax);
  
  if (ps_ffreal(*par,"cfl",&cfl)) 
    fprintf(stream,"NOTE: esg_readtimegrid - using default cfl = %e\n",cfl);
  
  if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
    fprintf(stream,"Error: esg_readtimegrid - either min dx=%e or cfl=%e "
	    " too small\n", a, cfl);
    return E_BADINPUT;
  }
  (model->tsind).dt = a*cfl/cmax;
  fprintf(stream, "esg_readtimegrid a = %g  cfl = %g  cmax = %g  dt = %g\n", a, cfl, cmax, (model->tsind).dt);
  return err;

}

static int esg_readmedia(PARARRAY * par, FILE * stream, IMODEL * model, int panelindex) {
  RDOM dom = model->ld_a;
  /* workspace for filename */
  char * pvelkey   = NULL;
  char * svelkey   = NULL;
  char * lambdakey = NULL;
  char * mukey     = NULL;
  char * rhokey    = NULL;
  char * denkey    = NULL;
  char * buoykey   = NULL;
  /**
   *  mp00 = lambda + 2 mu = vp * vp * rho
   *  mp01 = lambda = vp * vp * rho - 2 * vs * vs * rho
   *  ms0, ms1, ms2 = shifted mu
   *  rhoz, rhox, rhoy = shifted rho
   */
  int veloflag = 0; /*< flag that velocities was read */
  
  /* default values */
  ireal refpvel;
  ireal refsvel;
  ireal reflambda;
  ireal refmu;
  ireal refden;
  ireal refbuoy;

  /* parameter flags */
  int is_vel;    /* p-velocity and s-velocity is given */
  int is_lame;   /* Lame constants are given */
  int is_den;    /* density is given */
  int is_buoy;   /* buoyancy is given */
  
  char str_bn[20];   /* buoyancy */
  char str_dn[20];   /* density */
  char str_rh[20];   /* rho */
  
  char str_mu[20];   /* mu */
  
  char str_ind[2];   /* index of grid */

  /* tmp storage for read */
  RARR B0tmp, B1tmp;
  
  int i, j, n, m, l;
  int ndim, idim;
  IPNT ran;        /* length */
  IPNT rags;       /* start */
  IPNT rage;       /* end */
  IPNT rags_shifted; /* start for shifted array */
  IPNT rage_shifted;
  
  IPNT shift;
  int n_shift;
  int n_average;
  unsigned long ntot; /* total number of words in phys dom */
  ireal cmax = CMAX_DEF; /* max p-velocity from input */
  ireal cmin = CMIN_DEF; /* min p-velocity from input */
  ireal pvmax;            /* working pvmax */
  ireal pvmin;            /* working pvmin */

  int mult_mv;  /* multiple read flag for density or buoyancy */
  int mult_ms;  /* multiple read flag for mu or s-velocity */
    
  int err = 0;
  IPNT ii, kk;
  
  ireal q, b, q1, q2;
  IPNT gtype[RDOM_MAX_NARR];

  FD_MODEL * fdm = (FD_MODEL *)model->specs;
  

  /*--------------------------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/
  /* initialization */
  ra_setnull(&B0tmp);
  ra_setnull(&B1tmp);
  ndim = model->g.dim;
  
  esg_set_grid_type(stderr, ndim, gtype);
  /* must have a variable defined on PRIMAL grid in every axis */
  /* tedious, but works in general cases */
  err = rd_size(&dom, D_MP00, ran);
  if (err) {
    fprintf(stream, "Error: internal: allocated array [%s] size error #%d\n",
            fdm->ind2str(D_MP00), err);
    return E_INTERNAL;
  }
  err = rd_gse(&dom, D_MP00, rags, rage);
  if (err) {
    fprintf(stream, "Error: internal: allocated array [%s] gse error #%d\n",
            fdm->ind2str(D_MP00), err);
    return E_INTERNAL;
  }

  ntot = 1;
  for (idim = 0;idim < ndim;idim ++) ntot *= ran[idim];

  if (ps_ffreal(*par, "cmax", &cmax))
    fprintf(stream, "Note: esg_readmedia - using default cmax = %e\n", cmax);
  if (ps_ffreal(*par, "cmin", &cmin))
    fprintf(stream, "Note: esg_readmedia - using default cmin = %e\n", cmin);
  /*--------------------------------------------------------------------------*/
  /* reference values - first set defaults -----------------------------------*/
  /* lambda mu have the same unit with kappa ---------------------------------*/
  if ( strcmp(SEAM_UNITS, "SI") == 0 ) {
    reflambda = KAPPA_SI;
    refmu     = KAPPA_SI;
    refbuoy   = 1.0/RHO_SI;
  }
  else if ( strcmp(SEAM_UNITS, "MMSKG") == 0 ) {
    reflambda = KAPPA_MMSKG;
    refmu     = KAPPA_MMSKG;
    refbuoy   = 1.0/RHO_MMSKG;
  }
  else if ( strcmp(SEAM_UNITS, "KMSKG") == 0 ) {
    reflambda = KAPPA_KMSKG;
    refmu     = KAPPA_KMSKG;
    refbuoy   = 1.0/RHO_KMSKG;    
  }
  else if ( strcmp(SEAM_UNITS, "MMSGCM3") == 0 ) {
    reflambda = KAPPA_MMSGCM3;
    refmu     = KAPPA_MMSGCM3;
    refbuoy   = 1.0/RHO_MMSGCM3;
  }
  else {
    fprintf(stream, "Error: esg_readmedia from set default values\n");
    fprintf(stream, "SEAM units not set properly\n");
    return E_INTERNAL;
  }
  /*--------------------------------------------------------------------------*/
  /* detect which parameters are provided ------------------------------------*/
  refpvel = REAL_ZERO;
  refsvel = REAL_ZERO;
  refden  = REAL_ZERO;
    
  is_vel  = ( !ps_ffreal(*par, "refpvel", &refpvel) && refpvel > REAL_EPS );
  is_vel  = ( !ps_ffreal(*par, "refsvel", &refsvel) && refsvel > REAL_EPS ) && is_vel;
  is_lame = ( !ps_ffreal(*par, "reflambda", &reflambda) && reflambda > REAL_EPS );
  is_lame = ( !ps_ffreal(*par, "refmu", &refmu) && refmu > REAL_EPS ) && is_lame;
  is_den  = ( !ps_ffreal(*par, "refden", &refden) && refden > REAL_EPS );
  is_buoy = ( !ps_ffreal(*par, "refbuoy", &refbuoy) && refbuoy > REAL_EPS);
   
  if (!is_buoy) {
    if (is_vel && is_lame)  refbuoy = refsvel*refsvel/refmu;
    if (is_den)             refbuoy = REAL_ONE/refden;
  }
  if (!is_lame) {
    if (is_vel && is_den)  {
      refmu = refsvel*refsvel*refden;
      reflambda = refpvel*refpvel*refden - 2*refmu;
    }
    if (is_vel && is_buoy) {
      refmu = refsvel*refsvel/refbuoy;
      reflambda = refpvel*refpvel/refbuoy - 2*refmu;
    }
  }
    
  fprintf(stream, "Note: in esg_readmedia, reference values:\n");
  fprintf(stream, "  pvel = %e\n", refpvel);
  fprintf(stream, "  svel = %e\n", refsvel);
  fprintf(stream, "lambda = %e\n", reflambda);
  fprintf(stream, "    mu = %e\n", refmu);
  fprintf(stream, "   den = %e\n", refden);
  fprintf(stream, "  buoy = %e\n", refbuoy);

  /* check reference p-velocity for comformance */
  refpvel = sqrt((reflambda + 2*refmu)*refbuoy);
  if (refpvel < cmin || refpvel > cmax) {
    fprintf(stream, "Error: esg_readmedia - p-velocity out of bounds\n");
    return E_BADINPUT;
  }
  /**
   * File I/O section 
   * read Lame constants or p and s-velocity into D_MP00, D_MP01
   * since material parameters are static,  allocated domain is equal 
   * to computation domain, i.e., s0 == s
   */
  /* prerequisite: array mp00 and mp01 on primal grid in every direction */
  for(idim = 0;idim < ndim;idim ++) {
    if ( gtype[D_MP00][idim] != PRIMAL_GRID) {
      fprintf(stream, 
              "Error: esg_readmedia, array %s should be define on primal grid\n",
              fdm->ind2str(D_MP00));
      return E_BADINPUT;
    }
    if ( gtype[D_MP01][idim] != PRIMAL_GRID) {
      fprintf(stream, 
              "Error: esg_readmedia, array %s should be define on primal grid\n",
              fdm->ind2str(D_MP01));
      return E_BADINPUT;
    }
  }
  if ( !ps_ffcstring(*par, "lambda", &lambdakey) ) {
    if ( ps_ffcstring(*par, "mu", &mukey) ) {
      free(lambdakey);
      free(mukey);
      fprintf(stream, 
              "Error: in esg_readmedia: lambda and mu must be given at the same time\n");
      return E_BADINPUT;
    }
    /* lambda */
    err = rsfread(dom._s[D_MP01]._s0, rags, ran, lambdakey, 1, stream, panelindex);
    free(lambdakey);
    if (err) {
      fprintf(stream, "Error: esg_readmedia from rsfread - lambda\n");
      return err;
    }
    rd_gse(&dom, D_MP00, rags, rage);
    for (i=0; i<(ndim-1)*ndim/2; i++) {
      rd_gse(&dom,D_MS[i],rags_shifted,rage_shifted);
      for (j=0; j<ndim; j++) {
        rags[j] = iwave_min(rags_shifted[j], rags[j]);
        if (gtype[D_MS[i]][j] == DUAL_GRID) {
          rage[j] = iwave_max(rage_shifted[j] + 1, rage[j]);
        }
        else {
          rage[j] = iwave_max(rage_shifted[j], rage[j]);
        }
      }  
    }
    
    err = ra_create(&B1tmp, ndim, rags, rage);
    if (err) {
      fprintf(stream, "Error: esg_readmedia from rsfread - cannot allocate tmp0 array\n");
      return err;
    }
    ra_size(&B1tmp, ran);
    /* mu */
    err = rsfread(B1tmp._s, rags, ran, mukey, 1, stream, panelindex);
    if (err) {
      fprintf(stream,
              "Error: esg_readmedia from rsfread - mukey (single) = %s\n", mukey);
      ra_destroy(&B1tmp);
      return err;
    }
    free(mukey);
    
    rd_gse(&dom, D_MP00, rags, rage);
    for (idim = ndim;idim < RARR_MAX_NDIM; idim ++) {
      rags[idim] = 0;
      rage[idim] = 0;
    }
#if RARR_MAX_NDIM > 2
    for (ii[2] = rags[2];ii[2] <= rage[2];ii[2] ++) 
#endif
    {
#if RARR_MAX_NDIM > 1
      for (ii[1] = rags[1];ii[1] <= rage[1];ii[1] ++)
#endif
      {
        for (ii[0] = rags[0];ii[0] <= rage[0];ii[0] ++) {
          q = rd_gget(&dom, D_MP01, ii) + 2.0*ra_gget(&B1tmp, ii);
          rd_gset(&dom, D_MP00, ii, q);
        }
      }
    }
  }
  else if ( !ps_ffcstring(*par, "pvelocity", &pvelkey) ) {
    if ( ps_ffcstring(*par, "svelocity", &svelkey) ) {
      free(pvelkey);
      free(svelkey);
      fprintf(stream,
              "Error: in esg_readmedia: ");
      fprintf(stream, "pvelocity and svelocity must be given at the same time\n");
      return E_BADINPUT;
    }
    veloflag = 1;
    err = rsfread(dom._s[D_MP00]._s0, rags, ran, pvelkey, 1, stream, panelindex);
    free(pvelkey);
    if (err) {
      fprintf(stream, "Error: esg_readmedia from rsfread - pvelocity\n");
      return err;      
    }
    rd_gse(&dom, D_MP00, rags, rage);
    for (i=0; i<(ndim-1)*ndim/2; i++) {
      rd_gse(&dom,D_MS[i],rags_shifted,rage_shifted);
      for (j=0; j<ndim; j++) {
        rags[j] = iwave_min(rags_shifted[j], rags[j]);
        if (gtype[D_MS[i]][j] == DUAL_GRID) {
          rage[j] = iwave_max(rage_shifted[j] + 1, rage[j]);
        }
        else {
          rage[j] = iwave_max(rage_shifted[j], rage[j]);
        }
      }  
    }
    err = ra_create(&B1tmp, ndim, rags, rage);
    if (err) {
      fprintf(stream, "Error: esg_readmedia from rsfread - cannot allocate tmp0 array\n");
      return err;
    }
    ra_size(&B1tmp, ran);
    /* shear velocity */
    err = rsfread(B1tmp._s, rags, ran, svelkey, 1, stream, panelindex);
    if (err) {
      fprintf(stream,
              "Error: esg_readmedia from rsfread - svelkey (single) = %s\n", svelkey);
      ra_destroy(&B1tmp);
      return err;
    }
    free(svelkey);
    
    /* velocity check - performed only in case velocity is read directly */
    pvmax = dom._s[D_MP00]._s0[0];
    pvmin = pvmax;
    for (i = 0;i < ntot;i ++) {
      pvmax = iwave_max(pvmax, dom._s[D_MP00]._s0[i]);
      pvmin = iwave_min(pvmin, dom._s[D_MP00]._s0[i]);
    }
    err = err || (pvmin < cmin) || (pvmax > cmax);
    if (err) {
      fprintf(stream,
              "Error: esg_readmedia - pvelocity field read from ");
      fprintf(stream, "file violates bounds [%e, %e]\n", cmin, cmax);
      fprintf(stream, "min pvel = %e  max pvel = %e\n", pvmin, pvmax);
      rd_dump(&dom, D_MP00, stream);
      return E_OTHER;
    }
  }
  else {
    /* if no data read, set to reference values */
    for (i = 0;i < ntot;i ++) {
      dom._s[D_MP00]._s0[i] = reflambda + 2*refmu;
      dom._s[D_MP01]._s0[i] = reflambda;
    }
  }

  /* buoyancy - first look for multiple buoyancy */
  idim = 0;
  mult_mv = 1;
  while ( idim < ndim ) {
    strcpy(str_bn, "buoyancy");
    strcpy(str_dn, "density");
    strcpy(str_rh, "rho");
    sprintf(str_ind, "%d", idim + 1);
    strcat(str_bn, str_ind);
    strcat(str_dn, str_ind);
    strcat(str_rh, str_ind);

    if ( !ps_ffcstring(*par, str_bn, &buoykey) ||
         !ps_ffcstring(*par, str_dn, &denkey)  ||
         !ps_ffcstring(*par, str_rh, &rhokey) ) {
      if (veloflag) {
        fprintf(stream, "Error: esg_readmedia\n");
        fprintf(stream,
                "must provide Lame constants if combined with \n");
        fprintf(stream, "shifted density arrays\n");
        return E_FILE;
      }
    }
    else {
      mult_mv = 0;
      break;
    }
    rd_size(&dom, D_MV[idim], ran);
    rd_gse(&dom, D_MV[idim], rags, rage);
    if (buoykey) {
      err = rsfread(dom._s[D_MV[idim]]._s0, rags, ran, buoykey, 1, stream, panelindex);
      if (err) {
        fprintf(stream, 
                "Error: esg_readmedia from rsfread - buoykey = %s\n", buoykey);
        return err;
      }
    }
    else {
      if (denkey) {
        err = rsfread(dom._s[D_MV[idim]]._s0, rags, ran, denkey, 1, stream, panelindex);
        if (err) {
          fprintf(stream, 
                  "Error: esg_readmedia from rsfread - denkey = %s\n", denkey);
          return err;
        }
      }
      else if (rhokey) {
        err = rsfread(dom._s[D_MV[idim]]._s0, rags, ran, rhokey, 1, stream, panelindex);
        if (err) {
          fprintf(stream, 
                  "Error: esg_readmedia from rsfread - rhokey = %s\n", rhokey);
          return err;
        }
      }
      else {
        fprintf(stream, "Error: esg_readmedia from rsfread: ");
        fprintf(stream, "buoykey, denkey, rhokey are all empty\n");
        return E_BADINPUT;
      }
      ntot = 1;
      for (i = 0;i < ndim;i ++) ntot *= ran[i];
      for (i = 0;i < ntot;i ++)
        dom._s[D_MV[idim]]._s0[i] = REAL_ONE/(dom._s[D_MV[idim]]._s0[i]);      
    }
    if (buoykey)  {
      free(buoykey);
      buoykey = NULL;     
    }
    if (denkey) {
      free(denkey);
      denkey = NULL;
    }
    if (rhokey)  {
      free(rhokey);
      rhokey = NULL;
    }
    idim ++;
  }

  idim = 0;
  mult_ms = 1;
  mukey = NULL;
  while ( idim < ndim ) {
    strcpy(str_mu, "mu");
    sprintf(str_ind, "%d", idim + 1);
    strcat(str_mu, str_ind);

    if ( ps_ffcstring(*par, str_mu, &mukey) ) {
      mult_ms = 0;
      break;
    }
    rd_size(&dom, D_MS[idim], ran);
    rd_gse(&dom, D_MS[idim], rags, rage);
    ntot = 1;
    for (i = 0;i < ndim;i ++) ntot *= ran[i];
    err = rsfread(dom._s[D_MS[idim]]._s0, rags, ran, mukey, 1, stream, panelindex);
    fprintf(stream, "read mult mu%d: %s\n", idim+1, mukey);   
    free(mukey);
    mukey = NULL;
    if (err) {
      fprintf(stream, 
              "Error: esg_readmedia from rsfread - mukey = %s\n", mukey);
      return err;
    }
    idim ++;
  }

  /**
   * if fail to read multiple buoy or den names, read single field 
   * into workspace, then create shifted arrays by averaging
   */
  if (!mult_mv || veloflag) {
    buoykey = NULL;
    rhokey  = NULL;
    denkey  = NULL;
    if ( ps_ffcstring(*par, "buoyancy", &buoykey) &&
         ps_ffcstring(*par, "density", &denkey) &&
         ps_ffcstring(*par, "rho", &rhokey) ) {
      /* cannt read density or buoyancy, use ref value */
      fprintf(stream, "Warning: cann't read densit or buoyancy on primal grid, use ref value\n");
      for (i = 0;i < ndim;i ++) {
        rd_size(&dom, D_MV[i], ran);
        ntot = 1;
        for (j = 0;j < ndim;j ++) ntot *= ran[j];
        for (j = 0;j < ntot;j ++) dom._s[D_MV[i]]._s0[j] = refbuoy;
      }       
    }
    else {
      rd_gse(&dom, D_MP00, rags, rage);
      for (i=0; i<ndim; i++) {
        rd_gse(&dom,D_MV[i],rags_shifted,rage_shifted);
        for (j=0; j<ndim; j++) {
          rags[j] = iwave_min(rags_shifted[j], rags[j]);
          if (gtype[D_MV[i]][j] == DUAL_GRID) {
            rage[j] = iwave_max(rage_shifted[j] + 1, rage[j]);
          }
          else {
            rage[j] = iwave_max(rage_shifted[j], rage[j]);
          }
        }  
      }

      err = ra_create(&B0tmp, ndim, rags, rage);
      if (err) {
        fprintf(stream, "Error: esg_readmedia from rsfread - cannot allocate tmp0 array\n");
        return err;
      }
      ra_size(&B0tmp, ran);    
      if (buoykey) {
        err = rsfread(B0tmp._s, rags, ran, buoykey, 1, stream, panelindex);
        if (err) {
          fprintf(stream,
                  "Error: esg_readmedia from rsfread - buoykey (single) = %s\n", buoykey);
          ra_destroy(&B0tmp);
          return err;
        }
      }
      else {
        if (denkey) {
          err = rsfread(B0tmp._s, rags, ran, denkey, 1, stream, panelindex);
          if (err) {
            fprintf(stream, 
                    "Error: esg_readmedia from rsfread - denkey = %s\n", denkey);
            return err;
          }
        }
        else if (rhokey) {
          err = rsfread(B0tmp._s, rags, ran, rhokey, 1, stream, panelindex);
          if (err) {
            fprintf(stream, 
                    "Error: esg_readmedia from rsfread - rhokey = %s\n", rhokey);
            return err;
          }
        }
        else {
          fprintf(stream, "Error: esg_readmedia from rsfread: ");
          fprintf(stream, "buoykey, denkey, rhokey are all empty\n");
          return E_BADINPUT;
        }
        ntot = 1;
        for (i = 0;i < ndim;i ++) ntot *= ran[i];
        for (i = 0;i < ntot;i ++)
          B0tmp._s[i] = REAL_ONE/(B0tmp._s[i]);      
      }

      for (i = 0;i < ndim;i ++) {
        rd_gse(&dom, D_MV[i], rags, rage);
        for (j = ndim;j < RARR_MAX_NDIM;j ++) {
          rags[j] = 0;
          rage[j] = 0;
        }
        IASN(shift, IPNT_0);
        n_average = 1;
        n_shift = 0;
        for (j = 0;j < ndim;j ++) {
          if (gtype[D_MV[i]][j] == DUAL_GRID) {
            shift[j] = 1;
            n_average *= 2;
            n_shift ++;
          }
        }
        b = REAL_ONE/n_average;
#if RARR_MAX_NDIM > 2
        for (ii[2] = rags[2];ii[2] <= rage[2];ii[2] ++) 
#endif
        {
#if RARR_MAX_NDIM > 1
          for (ii[1] = rags[1];ii[1] <= rage[1];ii[1] ++)
#endif
          {
            for (ii[0] = rags[0];ii[0] <= rage[0];ii[0] ++) {
              //q = ra_gget(&B0tmp, ii);
              q = REAL_ONE / ra_gget(&B0tmp, ii);
              for (n = 1;n < n_average;n ++) {
                l = n;
                IASN(kk, ii);
                for (m = 0;m < ndim;m ++) {
                  if (shift[m]) {
                    kk[m] = ii[m] + l%2;
                    l = l/2 ;
                  }
                }
                //q += ra_gget(&B0tmp, kk);
                q += REAL_ONE / ra_gget(&B0tmp, kk);
              }
              //rd_gset(&dom, D_MV[i], ii, q*b);
              rd_gset(&dom, D_MV[i], ii, REAL_ONE/(q*b));
            }
          }
        }
      }
    }
    if (veloflag) {
      b = REAL_ONE/refbuoy;
      rd_gse(&dom, D_MP00, rags, rage);
      for (j = ndim;j < RARR_MAX_NDIM;j ++) {
        rags[j] = 0;
        rage[j] = 0;
      }
#if RARR_MAX_NDIM > 2
      for (ii[2] = rags[2];ii[2] <= rage[2];ii[2] ++) 
#endif
      {
#if RARR_MAX_NDIM > 1
        for (ii[1] = rags[1];ii[1] <=rage[1];ii[1] ++)
#endif
        {
          for (ii[0] = rags[0];ii[0] <= rage[0];ii[0] ++) {
            b = ra_gget(&B0tmp, ii); /* buoyancy */
            b = REAL_ONE/b;         /* density */
            q = rd_gget(&dom, D_MP00, ii);
            q = q*q*b;  /* lambda + 2 mu */
            rd_gset(&dom, D_MP00, ii, q);
            q1 = ra_gget(&B1tmp, ii);
            q1 = q1*q1*b; /* mu */
            rd_gset(&dom, D_MP01, ii, q-2*q1);            
          }
        }
      }
    }
    if (buoykey)  {
      free(buoykey);
      buoykey = NULL;     
    }
    if (denkey) {
      free(denkey);
      denkey = NULL;
    }
    if (rhokey)  {
      free(rhokey);
      rhokey = NULL;
    }
  }

  /**
   * if fail to read multiple mu, read single field 
   * into workspace, then create shifted arrays by averaging
   */
  if (!mult_ms) {
    mukey = NULL;
    svelkey = NULL;
    if ( ps_ffcstring(*par, "svelocity", &svelkey) &&
         ps_ffcstring(*par, "mu", &mukey) ) {
      /* cannt read mu or svelocity, use ref value */
      for (i = 0;i < ndim*(ndim-1)/2;i ++) {
        rd_size(&dom, D_MS[i], ran);
        ntot = 1;
        for (j = 0;j < ndim;j ++) ntot *= ran[j];
        for (j = 0;j < ntot;j ++) dom._s[D_MS[i]]._s0[j] = refmu;
      }
    }
    else {
      for (i = 0;i < ndim*(ndim-1)/2;i ++) {
        rd_gse(&dom, D_MS[i], rags, rage);
        for (j = ndim;j < RARR_MAX_NDIM;j ++) {
          rags[j] = 0;
          rage[j] = 0;
        }
        IASN(shift, IPNT_0);
        n_average = 1;
        n_shift = 0;
        for (j = 0;j < ndim;j ++) {
          if (gtype[D_MS[i]][j] == DUAL_GRID) {
            shift[j] = 1;
            n_average *= 2;
            n_shift ++;
          }
        }
        b = REAL_ONE/n_average;
#if RARR_MAX_NDIM > 2
        for (ii[2] = rags[2];ii[2] <= rage[2];ii[2] ++) 
#endif
        {
#if RARR_MAX_NDIM > 1
          for (ii[1] = rags[1];ii[1] <= rage[1];ii[1] ++)
#endif
          {
            for (ii[0] = rags[0];ii[0] <= rage[0];ii[0] ++) {
              if (veloflag) {
                q2 = ra_gget(&B1tmp, ii);  /* s-velocity */
                q1 = ra_gget(&B0tmp, ii);  /* buoyancy */
                q = q2*q2/q1;               /* mu */
              }
              else {
                q = ra_gget(&B1tmp, ii);  /* mu */
              }
              for (n = 1;n < n_average;n ++) {
                l = n;
                IASN(kk, ii);
                for (m = 0;m < ndim;m ++) {
                  if (shift[m]) {
                    kk[m] = ii[m] + l%2;
                    l = l/2;
                  }
                }
                if (veloflag) {
                  q2 = ra_gget(&B1tmp, kk);  /* s-velocity */
                  q1 = ra_gget(&B0tmp, kk); /* buoyancy */
                  q += q2*q2/q1;               /* mu */
                }
                else {
                  q += ra_gget(&B1tmp, kk);  /* mu */
                }
              }
              rd_gset(&dom, D_MS[i], ii, q*b);
            }
          }
        }
      }
    }
  }
  ra_destroy(&B0tmp);
  ra_destroy(&B1tmp);
  return 0;
}

/*----------------------------------------------------------------------------*/
/*   
     HERE I POPULATE ETA ARRAYS FROM THE MASTER INPUT FILE: 
*/
static int esgn_setetas(PARARRAY *pars, FILE *stream, IMODEL *model) {

  RDOM *dom;        /* allocated domain */
  int ndim;         /* number of dimensions */
  ireal pmlampl, L0, L1;                          /* pml variables */
  ireal cmax = CMAX_DEF;                          
  ireal dd, s;                                    /* dx/dy/dz */
  ireal eta0, eta1;                               /* pml etas */
  int idim, i, iarr;                             /* array index */
  IPNT gsa, gea, ind;
  ireal v0, v1;                                   /* physical bounds */
    
  dom = &(model->ld_a);
  ndim = model->g.dim;
    
  /* read cmax -------------------------------------------------------------*/
  if ( ps_ffreal(*pars, "cmax", &cmax) ) 
    fprintf(stream,"NOTE: esgn_setetas - using default cmax = %e\n",cmax);
  
  /* read pml amplitude ----------------------------------------------------*/
  pmlampl = 1.5 * log(1000.0);
  if ( ps_ffreal(*pars,"pml_ampl", &pmlampl) )
    fprintf(stream, "NOTE. Cannot read pml_ampl. Default = %g.\n",  pmlampl);
  else
    fprintf(stream, "NOTE. Eta multiplier pml_ampl = %g.\n", pmlampl);
  
  for ( idim = 0; idim < ndim; ++idim )
  {
    /* space step and area bounds */
    dd = model->g.axes[idim].d;
    v0 = 0.0;
    v1 = v0 + dd * (ireal)(model->g.axes[idim].n - 1);
    
    L0 = model->nls[idim] * dd; /* PML width in ireal length units */
    if ( L0 > 0.0 ) eta0 = pmlampl * cmax / L0;
    else            eta0 = 0.0;
    
    L1 = model->nrs[idim] * dd; /* PML width in ireal length units */
    if ( L1 > 0.0 ) eta1 = pmlampl * cmax / L1;
    else            eta1 = 0.0;
    
    iarr = D_EP[idim];
    rd_gse(dom, iarr, gsa, gea);
    for ( i = 1; i < RARR_MAX_NDIM; ++i ) ind[i] = gsa[i] = gea[i] = 0;
    
    for ( ind[0] = gsa[0]; ind[0] <= gea[0]; ++ind[0] ) 
    {
      s = v0 + dd * (ireal)(ind[0]);
      if      ( s < v0 ) 
      {
        s = (s - v0) / L0;
        s = eta0 * fabs(s * s);
      }
      else if ( s > v1 )
      {
        s = (s - v1) / L1;
        s = eta1 * fabs(s * s);
      }
      else
      {
        s = 0.0;
      }
      rd_gset(dom, iarr, ind, s);
    }
    
    iarr = D_EV[idim];
    rd_gse(dom, iarr, gsa, gea);
    for ( i = 1; i < RARR_MAX_NDIM; ++i ) ind[i] = gsa[i] = gea[i] = 0;
    
    for ( ind[0] = gsa[0]; ind[0] <= gea[0]; ++ind[0] ) 
    {
      s = v0 + dd * (0.5 + (ireal)(ind[0]));
      if      ( s < v0 ) 
      {
        s = (s - v0) / L0;
        s = eta0 * fabs(s * s);
      }
      else if ( s > v1 )
      {
        s = (s - v1) / L1;
        s = eta1 * fabs(s * s);
      }
      else
      {
        s = 0.0;
      }
      rd_gset(dom, iarr, ind, s);
    }        
  }
  return 0;
}    

/* new interface for readmedia */
int esgn_readmedia(PARARRAY * pars, 
		  FILE * stream,
		  IMODEL * model,
		  int panelindex) {

  int err=0;

  /* read main medium files */
  err=esg_readmedia(pars,
		   stream,
		   model,
		   panelindex);

  if (!err) err=esgn_setetas(pars,stream,model);
  return err;
}

/* helper function - not in header file */
int esg_readpmlgrid(IPNT nl, IPNT nr, 
		    RPNT dx, size_t ndim, 
		    PARARRAY par, FILE * stream) {
  int idim;
  ireal fpeak = FPEAK_DEF;
  ireal cmax  = CMAX_DEF;
  ireal tmp;
  char key[30];
  size_t kl=4;

  if (ps_ffreal(par,"fpeak",&fpeak))
    fprintf(stream,"NOTE: sgn_readpmlgrid - using default fpeak = %e\n",fpeak);

  if (ps_ffreal(par,"cmax",&cmax)) 
    fprintf(stream,"NOTE: sgn_readpmlgrid - using default cmax = %e\n",cmax);

  for (idim=0; idim<ndim; idim++ ) {
	  if (dx[idim]<REAL_EPS) {
		  fprintf(stream, "dx[%d] = %f\n", idim, dx[idim]);
		  return E_BADINPUT;
	  }
    snprintf(key,kl,"nl%d",idim+1);
    tmp=0.0;
    ps_ffreal(par, key, &tmp);

    nl[idim]=iwave_max(0,(int)((ceil)(tmp*cmax/(fpeak*dx[idim]))));
    snprintf(key,kl,"nr%d",idim+1);
    tmp=0.0;
    ps_ffreal(par,key,&tmp);
    nr[idim]=iwave_max(0,(int)((ceil)(tmp*cmax/(fpeak*dx[idim]))));
  }
  return 0;
}

/*----------------------------------------------------------------------------*/
/* 
  Number of pressure grid points in physical domain and PML domains.
  sgnm->ns  has IPNT type and contains physical domain size.
  sgnm->nsl has IPNT type and contains left PML domains' sizes.
  sgnm->nsr has IPNT type and contains right PML domains' sizes.
  dxs has RPNT type and contains space steps.
  each element (0,1,2) of IPNT/RPNT corresponds to a space dimension.

  Signature change 10.07.08 WWS: removed dt, not set here 

  Major sig change 10.11.10 WWS: adapt to XW's readgrid interface

*/
int esg_readgrid(PARARRAY * pars, 
		 FILE *stream, 
		 IMODEL * model) {

  // workspace
  int err=0;
  RPNT d;
	
  IASN(model->nls, IPNT_0);  /* default left PML width */
  IASN(model->nrs, IPNT_0);  /* default right PML width */
  
  /* read physical grid parameters from bulk modulus or velocity
     grid, or from parameter file */
  if ((err=init_elastic_geom_par(&(model->g),*pars,stream))) {
      fprintf(stream,"Error: could not read grid parameters\n");
      return err;
  }

  /* extract grid steps */
  get_d(d,model->g);

  /* extract pml sizes from parameter files */
  err=esg_readpmlgrid(model->nls,model->nrs,d, (model->g).dim, *pars,stream);
  if (err) {
    fprintf(stream,"Error: could not read pml grid parameters\n");
    return err;
  }
  
  return err;
}

