/* interface to sg_readtimegrid */

#include "fd.h"
#include "sgnpars.h"
#include "sgn_read.h"

int set_default_material_params(ireal * refkappa,
				ireal * refbuoy) {
  if ( strcmp(SEAM_UNITS, "SI") == 0 ) {
    *refkappa=KAPPA_SI;
    *refbuoy =1.0/RHO_SI;
  }
  else if ( strcmp(SEAM_UNITS, "MMSKG") == 0 ) {
    *refkappa=KAPPA_MMSKG;
    *refbuoy =1.0/RHO_MMSKG;
  }
  else if ( strcmp(SEAM_UNITS, "KMSKG") == 0 ) {
    *refkappa=KAPPA_KMSKG;
    *refbuoy =1.0/RHO_KMSKG;    
  }
  else if ( strcmp(SEAM_UNITS, "MMSGCM3") == 0 ) {
    *refkappa=KAPPA_MMSGCM3;
    *refbuoy =1.0/RHO_MMSGCM3;
  }
  else {
    return 1;
  }
  return 0;
}

/* this version: reads grid geometry from bulk modulus or velocity RSF
   file, or from parameter list. */

int init_acoustic_geom_par(grid * _g, PARARRAY par, FILE * fp) {

  char * _kappakey=NULL;
  char * _bmodkey=NULL;
  char * _velokey =NULL;
  int err=0;

  int i;
  int tmp;
  char key[3];
  size_t kl=3;

  /* check for permissible velo/bulkmod keywords */  
  ps_flcstring(par,"kappafile",&_kappakey);
  ps_flcstring(par,"bulkmod",&_bmodkey);
  ps_flcstring(par,"velocity",&_velokey);

  if (_kappakey) {
    err=read_grid(_g,_kappakey,fp);
    if (_kappakey) userfree_(_kappakey);
    if (err) return err;
  }
  else if (_bmodkey) {
    err=read_grid(_g,_bmodkey,fp);
    if (_bmodkey) userfree_(_bmodkey);
    if (err) return err;
  }
  else if (_velokey) {
    err=read_grid(_g,_velokey,fp);
    if (_velokey) userfree_(_velokey);
    if (err) return err;
  }
  else {
    /* look for n's, d's, and o's in parameter list */
    _g->dim=0;
    for (i=0;i<RARR_MAX_NDIM;i++) {
      snprintf(key,kl,"n%d",i+1);
      tmp=1;
      ps_flint(par,key,&tmp);
      _g->axes[i].n=tmp;
      snprintf(key,kl,"d%d",i+1);
      if (ps_flreal(par,key,&(_g->axes[i].d))) _g->axes[i].d=1.0;
      snprintf(key,kl,"o%d",i+1);
      if (ps_flreal(par,key,&(_g->axes[i].o))) _g->axes[i].o=0.0;
      /* determine dim by finding least axis index with n>1 */
      if (_g->axes[i].n>1) _g->dim=iwave_max(_g->dim,i+1);
    }

    /* added 03.03.09 - need to determine order params! */
    /*    if (_g->dim > 0) { */
    if (RARR_MAX_NDIM > 0) {
      tmp=1;
      ps_flint(par,"z_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>RARR_MAX_NDIM-1) {
	fprintf(fp,"ERROR: init_acoustic_geom\n");
	fprintf(fp,"z_axis index = %d out of range for dim = %d\n",
		tmp,RARR_MAX_NDIM);
	return E_OTHER;
      }
      _g->axes[tmp].id=0;
    }
    /*    if (_g->dim > 1) { */
    if (RARR_MAX_NDIM > 1) {
      tmp=2;
      ps_flint(par,"x_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>RARR_MAX_NDIM-1) {
	fprintf(fp,"ERROR: init_acoustic_geom\n");
	fprintf(fp,"x_axis index = %d out of range for dim = %d\n",
		tmp,RARR_MAX_NDIM);
	return E_OTHER;
      }
      _g->axes[tmp].id=1;
    }

    /*    if (_g->dim > 2) {*/
    if (RARR_MAX_NDIM > 2) {
      tmp=3;
      ps_flint(par,"y_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>RARR_MAX_NDIM-1) {
	fprintf(fp,"ERROR: init_acoustic_geom\n");
	fprintf(fp,"y_axis index = %d out of range for dim = %d\n",
		tmp,RARR_MAX_NDIM);
	return E_OTHER;
      }
      _g->axes[tmp].id=2;
    }
  }

  return err;
}

/* D.S., Feb 19, 2010. compute dt only based on cfl value and maximum
   velocity provided in para file */
int sg_readcfltime(RDOM dom, 
		   FILE *stream, 
		   ireal *dt, 
		   grid g, 
		   ireal cfl, 
		   ireal cmax,
		   PARARRAY par) {
  
  int i;                        /* counter */
  int ndim;                     /* model dimn */
  int err = 0;                  /* error flag */
  ireal a;              
  /*
  char * mode = NULL;
  */
  
  ndim = g.dim;
  
  /*
  ps_flcstring(par,"mode",&mode);      
  */
  a = g.axes[0].d;
  
  /*
  if (mode==NULL)
  */
  for ( i = 1; i < ndim; ++i ) a = iwave_min(a, g.axes[i].d);
  
  if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
    fprintf(stream,"ERROR: sg_readcfltime - either min dx=%e or cfl=%e "
	    " too small\n", a, cfl);
    return E_BADINPUT;
  }
  *dt = a*cfl/(cmax*sqrt((float)(ndim)));
  
  return err;

}

/* computes stable dt according to Gustafsson and Wahlund, SIAM
   J. Sci. Comput. vol. 26, 2004, p.264 generalized to 3D.
 */
 
int sg_readgustime(RDOM dom, 
		   FILE *stream, 
		   ireal *dt, 
		   grid g, 
		   float cfl, 
		   PARARRAY par) {
  IPNT kgs, kge, bgs, bge, ip;
  int err, ndim, i, j;
  ireal cent, right, left, vgus, a, tmp;
  ireal cp[3];
  ireal cv[3];
#ifdef IWAVE_USE_MPI
  MPI_Comm cm;
#endif
  /* WWS disabled 02.10
    char * mode = NULL; */                       /* DS 04.05.09 */
  
  err = 0;
  
  ndim = g.dim;
  /* set extra dimensions to 0 */
  for ( i = ndim; i < RARR_MAX_NDIM; ++i ){
    kgs[i] = kge[i] = 0;	
    bgs[i] = bge[i] = 0;	
  }
  
  rd_gse(&dom, D_MP0, kgs, kge);
  /* compute cp[i] = max_{i,j,k} 0.5 * ( sqrt(kappa_i * buoyancy_i-1) + 
     sqrt(kappa_i * buoyancy_i) )
  */
  for (i = 0; i < ndim; ++i){
    rd_gse(&dom, D_MV[i], bgs, bge);
    
    /* push bulk modulus bdr inside buoyancy bdr */
    for (j = 0; j < ndim; ++j)
      if (i == j){
	if ( kgs[j] <= bgs[j] ) kgs[j] = bgs[j] + 1;
	if ( kge[i] >  bge[j] ) kge[j] = bge[j];
      }
      else {
	if ( kgs[j] <  bgs[j] ) kgs[j] = bgs[j];
	if ( kge[i] >  bge[j] ) kge[j] = bge[j];
	
      }
    
    cp[i] = 0.0;
#if RARR_MAX_NDIM > 2
    for ( ip[2] = kgs[2]; ip[2] <= kge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
      for ( ip[1] = kgs[1]; ip[1] <= kge[1]; ++ip[1] )
#endif
	for ( ip[0] = kgs[0]; ip[0] <= kge[0]; ++ip[0] ){
	  
	  cent = rd_gget(&dom, D_MP0, ip);
	  right = rd_gget(&dom, D_MV[i], ip);
	  ip[i]--;
	  left = rd_gget(&dom, D_MV[i], ip);
	  ip[i]++;
	  cp[i] = iwave_max(cp[i], 0.5 * (sqrt(cent * right) + sqrt(cent * left)) );
	}	
  }
  
  /* bulk modulus bdr might be corrupted, read again */
  rd_gse(&dom, D_MP0, kgs, kge);
  /* compute cv[i] = max_{i,j,k} 0.5 * ( sqrt(kappa_i * buoyancy_i-1) + 
     sqrt(kappa_i-1 * buoyancy_i-1) )
  */
  for (i = 0; i < ndim; ++i){
    rd_gse(&dom, D_MV[i], bgs, bge);
    
    for (j = 0; j < ndim; ++j)
      if (i == j){
	if ( bgs[j] <  kgs[j] ) bgs[j] = kgs[j];
	if ( bge[j] >= kge[j] ) bge[j] = kge[j] - 1;
      }
      else {
	if ( bgs[j] <  kgs[j] ) bgs[j] = kgs[j];
	if ( bge[j] >  kge[j])  bge[j] = kge[j];
	
      }
    
    cv[i] = 0.0;
#if RARR_MAX_NDIM > 2
    for ( ip[2] = bgs[2]; ip[2] <= bge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
      for ( ip[1] = bgs[1]; ip[1] <= bge[1]; ++ip[1] )
#endif
	for ( ip[0] = bgs[0]; ip[0] <= bge[0]; ++ip[0] ){
	  
	  cent = rd_gget(&dom, D_MV[i], ip);
	  left = rd_gget(&dom, D_MP0, ip);
	  ip[i]++;
	  right = rd_gget(&dom, D_MP0, ip);
	  ip[i]--;
	  cv[i] = iwave_max(cp[i], 0.5 * (sqrt(right * cent) + sqrt(left * cent)) );
	}	
    /*fprintf(stderr, "cv[%d] = %g\n\n\n", i, cv[i]);*/
  }
  
  tmp = vgus =  0.0;
  /* DS 04.05.09 */
  /* disabled - WWS 02.10 */
  /*
  ps_flcstring(par,"mode",&mode);      
  if (mode)
	tmp = iwave_max( tmp, cp[0] );
  else
  */
  for (i = 0; i < ndim; i++) tmp = iwave_max( tmp, sqrt(cp[i] * cv[i]) );
  
#ifdef IWAVE_USE_MPI
  cm = retrieveComm();
  MPI_Allreduce(&tmp, &vgus, 1, IWAVE_MPI_REAL, MPI_MAX, cm);
#else
  vgus = tmp;
#endif

  /* get d, compute dt */
  a = g.axes[0].d;	
  /*  if (mode==NULL) */   /* DS 04.05.09 */
  for ( i = 1; i < ndim; ++i ) a = iwave_min(a, g.axes[i].d);
  
  if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
    fprintf(stream,"ERROR: sg_readgustime - either min dx=%e or cfl=%e "
	    " too small\n", a, cfl);
    return E_BADINPUT;
  }
  
  *dt = a*cfl/(vgus*sqrt((float)(ndim)));
  fprintf(stream, "NOTE - sg_readgustime: Computed dt=%12.4e from cfl, bulk modulus, density fields\n", *dt);
  fprintf(stream, "NOTE. min space step=%e, cfl=%e, Gustafsson factor=%e\n",a,cfl,vgus);

  return err;
}

/* implements new time grid logic:
   - by default choose stable time step (max_step set), optionally with 
     reduced cfl
   - if max_step unset, use CFL timestep but test against max step
   WWS 02.10
   NEW INTERFACE 22.11.10 WWS
   old interface
   int sg_readtimegrid(RDOM dom, FILE *stream, ireal *dt, grid g, PARARRAY par) {
   where
   dom = model->ld_a;
   dt  = &((model->tsind).dt);
   g=model->g;
   par=*pars;

*/

int asg_readtimegrid(PARARRAY *pars, FILE * stream, IMODEL * model) {

  ireal cmax = CMAX_DEF;        /* from input */
  ireal cfl = CFL_DEF;          /* default cfl fraction */
  ireal cflgus = CFL_DEF;       /* cfl frac to compute max stable step */
  ireal dtgus;                  /* max stable dt */
  int max_step = 1;             /* flag to take max permitted step */
  int err=0;
  RDOM dom;                     /* allocated domain */
  ireal * dt;                   /* time step, to be computed */
  grid g;                       /* spatial grid */
  PARARRAY par;                 /* parameter array */

  /* assign variables from old interface */
  dom = model->ld_a;
  dt  = &((model->tsind).dt);
  g=model->g;
  par=*pars;

  /* branch on presence of parameter dt */
  if ( !ps_flreal(par,"dt", dt ) ){
    fprintf(stream,"NOTE: sg_readtimegrid - dt=%12.4e read from param table\n", *dt);	
    fprintf(stream,"NOTE: NOT CHECKED FOR STABILITY!\n");
    return 0;
  }

  if (ps_flreal(par,"cmax",&cmax)) 
    fprintf(stream,"NOTE: sg_readtimegrid - using default cmax = %e\n",cmax);
  
  if (ps_flreal(par,"cfl",&cfl)) {
    fprintf(stream,"NOTE: sg_readtimegrid - using default cfl = %e\n",cfl);
  }
  else {
    if (cfl>CFL_DEF || cfl<REAL_EPS) {
      fprintf(stream,"ERROR: sg_readtimegrid - input cfl fraction = %e exceeds default = %e or too small\n",cfl,CFL_DEF);
      return E_BADINPUT;
    }
  }
    
	
  ps_flint(par,"max_step",&max_step);		
  if (max_step) {
    fprintf(stream,"NOTE: sg_readtimegrid - dt computed from max stable step, CFL fraction = %e\n",cfl);
    cflgus=cfl;
  }
  else {
    fprintf(stream,"NOTE: sg_readtimegrid - base dt on asserted CMAX = %e, CFL fraction = %e\n",cmax,cfl);
  }

  /* compute max stable step, optionally scaled by cfl from table */
  if (err=sg_readgustime(dom,stream,&dtgus,g,cflgus,par)) {  
    fprintf(stream,"\ncalled from sg_readtimegrid\n");
    return err;
  }

  if (max_step) {
    *dt=dtgus;
    return 0;
  }

  if (err=sg_readcfltime(dom,stream,dt,g,cfl,cmax,par)) {  
    fprintf(stream,"\ncalled from sg_readtimegrid\n");
    return err;
  }
  
  if (*dt > dtgus) {
    fprintf(stream,"ERROR: sg_readtimegrid - max_step option\n");
    fprintf(stream,"computed dt based on CFL criterion only = %e\n",*dt);
    fprintf(stream,"exceeds max stable step = %e\n",dtgus);
    return E_BADINPUT;
  }

  return 0;
}

/** Input mechanical parameters. Target: bulk modulus and buoyancy on
    appropriately shifted grids. Multiple possible sources: bulk
    modulus or velocity on integer grid, multiple density or buoyancy
    fields on half-integer shifted grids, single density or buoyancy
    on integer grids. Both classes of parameters defaulted to
    appropriate homogeneous medium values for water. 

    Includes velocity check.

    * @param[in]   panelindex (int) - panel index of extended model (always be 0 for non-extended model) 
    
    22.11.10: WWS
    This is the old interface from the pre-version 2 IWAVE.
    couple to function which sets eta variables to implement
    new interface
*/

int sg_readmedia(RDOM dom,
		 FILE * stream,
		 PARARRAY par,
		 int panelindex) {

  /*==================================================
    ================= DECLARATIONS ===================
    ==================================================*/

  /* workspace for filenames */  
  char * kappakey=NULL;
  char * bulkkey=NULL;
  char * velokey=NULL;
  int veloflag = 0; /* flag that velo was read */
  char * rhokey=NULL;
  char * buoykey=NULL;

  /* default values */
  ireal refkappa;
  ireal refbuoy;
  ireal refvel;
  ireal refden;

  /* parameter flags */
  int is_kap;
  int is_bou;
  int is_vel;
  int is_den;

  /* workspace */
  char bn[20];
  char dn[20];
  char rh[20];
  char dg[2];

  /* tmp storage for read */
  RARR Btmp;

  /* other workspace */
  int i,j,n,m; /* counters */
  int dim;         /* problem dimension */
  IPNT ran;        /* axis lens */
  IPNT rags;       /* axis starts */
  IPNT rage;       /* axis ends */
  unsigned long ntot;        /* total num of words in phys dom */
  ireal cmax = CMAX_DEF;  /* max velo */
  ireal cmin = CMIN_DEF;  /* min velo */
  ireal dmax = DMAX_DEF;  /* max dens */
  ireal dmin = DMIN_DEF;  /* min dens */
  /* bulk modulus, buoyancy limits computed from velo, density */
  ireal kmax;
  ireal kmin;
  ireal bmax;
  ireal bmin;

  ireal vmax;  /* working vmax */
  ireal vmin;  /* working vmin */

  int mult;   /* flags multiple dn/buoy input */
  /* error flag */
  int err=0;
  IPNT ip;
  IPNT raVgs, raVge, raPgs, raPge;
  ireal q, b;
  
  /*==================================================
    =============== END DECLARATIONS =================
    ==================================================*/  

  /* initializations */

  ra_setnull(&Btmp);
  /* axis lens, starts for bulk mod - allocated, not virtual */
  rd_ndim(&dom,D_MP0,&dim);
  rd_size(&dom,D_MP0,ran);
  rd_gse(&dom, D_MP0,rags,rage);
  ntot=1;
  for (i=0;i<dim;i++) ntot*=ran[i];

  /* SANITY CHECK MOVED TO SINGLE BOYANCY - I&T 03/05 */
  
  if (ps_flreal(par,"cmax",&cmax)) 
    fprintf(stream,"NOTE: sg_readmedia - using default max velocity = %e\n",cmax);

  if (ps_flreal(par,"cmin",&cmin)) 
    fprintf(stream,"NOTE: sg_readmedia - using default min velocity = %e\n",cmin);

  if (ps_flreal(par,"dmax",&dmax)) 
    fprintf(stream,"NOTE: sg_readmedia - using default max density = %e\n",dmax);

  if (ps_flreal(par,"dmin",&dmin)) 
    fprintf(stream,"NOTE: sg_readmedia - using default min density = %e\n",dmin);

  /* set limits for bulk mod, buoyancy */
  kmax=cmax*cmax*dmax;
  kmin=cmin*cmin*dmin;
  bmax=REAL_ONE/dmin;
  bmin=REAL_ONE/dmax;

  /* reference values - first set defaults */
  if (err=set_default_material_params(&refkappa,
				      &refbuoy)) {
    fprintf(stream,"ERROR: sg_readmedia from set_default\n");
    fprintf(stream,"SEAM units not set properly\n");
    return err;
  }

  /* WWS 05.12.08 - then look for reference values of any valid combination
     of parameters in the param list 
  */
  /* detect which parameters are provided */
  refvel=REAL_ZERO;
  refden=REAL_ZERO;
  is_kap=(!ps_flreal(par,"refkappa",&refkappa) &&(refkappa > REAL_EPS));
  is_bou=(!ps_flreal(par,"refbuoy",&refbuoy) && (refbuoy > REAL_EPS));
  is_vel=(!ps_flreal(par,"refvel",&refvel) && (refvel > REAL_EPS));
  is_den=(!ps_flreal(par,"refden",&refden) && (refden > REAL_EPS));
  /* if buoyancy is not in param list, but either density or velocity 
     are, then compute - otherwise use default value */
  if (!is_bou) {
    /* buoyancy computations */
    if (is_kap && is_vel) refbuoy=refvel*refvel/refkappa;
    if (is_den) refbuoy=REAL_ONE/refden;
  }
  /* if bulk mod is not in param list, but (density or buoyancy) and 
     velocity are, then comput - otherwise use default value */
  if (!is_kap) {
    if (is_vel && is_den) refkappa=refden*refvel*refvel;
    if (is_vel && is_bou) refkappa=refvel*refvel/refbuoy;
  }
  fprintf(stream,
	  "NOTE: in sg_readmedia, reference values: \n");
  fprintf(stream,
	  "  kap = %e\n",refkappa);
  fprintf(stream,
	  "  bou = %e\n",refbuoy);
  fprintf(stream,
	  "  vel = %e\n",refvel);
  fprintf(stream,
	  "  den = %e\n",refden);
  /* at this point reference values of bulk mod and buoyancy should 
     be available. check reference velocity for conformance */
  refvel=sqrt(refkappa*refbuoy);
  if ((refvel<cmin) || (refvel>cmax)) {
    fprintf(stream,
	    "ERROR: sg_readmedia - reference velocity out of bounds\n");
    return E_BADINPUT;
  } 

  /* FILE I/O SECTION */

  /* read bulkmod or velocity into D_MP0 - flag velocity */
  if (!ps_flcstring(par,"kappa", &kappakey) ) {
    if (!ps_flcstring(par,"bulkmod", &bulkkey)) {
      userfree_(kappakey);
      userfree_(bulkkey);
      fprintf(stream,
	      "ERROR: sg_readmedia from ra_rsfread: both kappa and bulkmod"
	      " are supplied (ambiguity)\n");
      return E_BADINPUT;
    }
		  
    err=rsfread(dom._s[D_MP0]._s0,rags,ran,kappakey,1,stream, panelindex);
    userfree_(kappakey);
    if (err) {
      fprintf(stream,"ERROR: sg_readmedia from ra_rsfread - kappa\n");
      return err;
    }
    /* bulk mod check */
    vmax=dom._s[D_MP0]._s0[0];
    vmin=dom._s[D_MP0]._s0[0];
    for (j=1;j<ntot;j++) {
      vmax=iwave_max(vmax,dom._s[D_MP0]._s0[j]);
      vmin=iwave_min(vmin,dom._s[D_MP0]._s0[j]);
    }
    err = (vmin < kmin) || (vmax > kmax);
    if (err) {
      fprintf(stream,
	      "ERROR: sg_readmedia - bulk mod field read from file violates bounds [%e, %e]\n",
	      kmin,kmax);
      fprintf(stream,"min kappa = %e max kappa = %e\n",vmin,vmax);
      return E_OTHER;
    }
  }
  else if (!ps_flcstring(par,"bulkmod", &bulkkey)) {
    err=rsfread(dom._s[D_MP0]._s0,rags,ran,bulkkey,1,stream,panelindex);
    userfree_(bulkkey);
    if (err) {
      fprintf(stream,"ERROR: sg_readmedia from ra_rsfread - bulkmod\n");
      return err;
    }
    /* bulk mod check */
    vmax=dom._s[D_MP0]._s0[0];
    vmin=dom._s[D_MP0]._s0[0];
    for (j=1;j<ntot;j++) {
      vmax=iwave_max(vmax,dom._s[D_MP0]._s0[j]);
      vmin=iwave_min(vmin,dom._s[D_MP0]._s0[j]);
    }
    err = (vmin < kmin) || (vmax > kmax);
    if (err) {
      fprintf(stream,
	      "ERROR: sg_readmedia - bulk mod field read from file violates bounds [%e, %e]\n",
	      kmin,kmax);
      fprintf(stream,"min kappa = %e max kappa = %e\n",vmin,vmax);
      return E_OTHER;
    }
  }
  else if (!ps_flcstring(par,"velocity",&velokey)) {
    fprintf(stream,"sg_readmedia -> rsfread velo\n");
    fflush(stream);
    err=rsfread(dom._s[D_MP0]._s0,rags,ran,velokey,1,stream,panelindex);
    fprintf(stream,"sg_readmedia <- rsfread velo\n");
    fflush(stream);
    userfree_(velokey);
    veloflag = 1;
    if (err) {
      fprintf(stream,"ERROR: sg_readmedia from ra_rsfread - velocity\n");
      return err;
    }

    /* velocity check - performed only in case velocity is read directly */
    vmax=dom._s[D_MP0]._s0[0];
    vmin=dom._s[D_MP0]._s0[0];
    for (j=1;j<ntot;j++) {
      vmax=iwave_max(vmax,dom._s[D_MP0]._s0[j]);
      vmin=iwave_min(vmin,dom._s[D_MP0]._s0[j]);
    }
    err = (vmin < cmin) || (vmax > cmax);
    if (err) {
      fprintf(stream,
	      "ERROR: sg_readmedia - velocity field read from file violates bounds [%e, %e]\n",
	      cmin,cmax);
      fprintf(stream,"min v = %e max v = %e\n",vmin,vmax);
      return E_OTHER;
    }

  }
  else {
    /* if no data read, set to ref value */
    for (j=0;j<ntot;j++) dom._s[D_MP0]._s0[j]=refkappa;
  }

  /* buoyancy - first look for multiple buoyancy */
  i=0;
  mult=1;
  while ((i<dim) && mult) {
    strcpy(bn,"buoyancy");
    strcpy(dn,"density");
    strcpy(rh,"rho");
    sprintf(dg,"%d",i + 1);
    strcat(bn,dg);
    strcat(dn,dg);
    strcat(rh,dg);

    if (!ps_flcstring(par,bn,&buoykey)) {
      if (veloflag) {
	fprintf(stream,"ERROR: sg_readmedia\n");
	fprintf(stream,
		"must provide bulk modulus directly if combined with \n");
	fprintf(stream,"shifted buoyancy arrays\n");
	return E_FILE;
      }
      /* modify axis lengths for each buoyancy */
      rd_size(&dom,D_MV[i],ran);
      rd_gse(&dom,D_MV[i],rags,rage);
      err = rsfread(dom._s[D_MV[i]]._s0,rags,ran,buoykey,1,stream,panelindex);
      userfree_(buoykey);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia from ra_rsfread - buoykey=%s\n",buoykey);
	return err;
      }
      /* buoyancy check */
      vmax=dom._s[D_MV[i]]._s0[0];
      vmin=dom._s[D_MV[i]]._s0[0];
      for (j=1;j<ntot;j++) {
	vmax=iwave_max(vmax,dom._s[D_MV[i]]._s0[j]);
	vmin=iwave_min(vmin,dom._s[D_MV[i]]._s0[j]);
      }
      err = (vmin < bmin) || (vmax > bmax);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia - buoyancy field %d read from file violates bounds [%e, %e]\n",i,bmin,bmax);
	fprintf(stream,"min buoy = %e max buoy = %e\n",vmin,vmax);
	return E_OTHER;
      }

    }
    else if ((!ps_flcstring(par,dn,&rhokey)) || 
	     (!ps_flcstring(par,rh,&rhokey))) {  /* TODO: Memory leak in rhokey here !!!!!! */
      if (veloflag) {
	fprintf(stream,"ERROR: sg_readmedia\n");
	fprintf(stream,
		"must provide bulk modulus directly if combined with \n");
	fprintf(stream,"shifted density arrays\n");
	return E_FILE;
      }
      rd_size(&dom,D_MV[i],ran);
      rd_gse(&dom,D_MV[i],rags,rage);
      err=rsfread(dom._s[D_MV[i]]._s0,rags,ran,rhokey,1,stream,panelindex);
      userfree_(rhokey);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia from ra_rsfread - rhokey=%s\n",rhokey);
	return err;
      }              
      /* density check */
      vmax=dom._s[D_MV[i]]._s0[0];
      vmin=dom._s[D_MV[i]]._s0[0];
      for (j=1;j<ntot;j++) {
	vmax=iwave_max(vmax,dom._s[D_MV[i]]._s0[j]);
	vmin=iwave_min(vmin,dom._s[D_MV[i]]._s0[j]);
      }
      err = (vmin < dmin) || (vmax > dmax);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia - desnsity field %d read from file violates bounds [%e, %e]\n",i,dmin,dmax);
	fprintf(stream,"min density = %e max density = %e\n",vmin,vmax);
	return E_OTHER;
      }

      n=1;
      for (j=0;j<dim;j++) n*=ran[j];
      /* THIS SHOULD BE PROTECTED SOMEHOW AGAINS ZERODIVIDE */
      for (j=0;j<n;j++) dom._s[D_MV[i]]._s0[j] 
	= REAL_ONE/(dom._s[D_MV[i]]._s0[j]);
    }
    else {
      mult=0;
    }
    i++;
  }
  
  /* if we failed to read multiple buoy or den names, read single
     field into workspace, create shifted arrays by averaging. */

  if (!mult) {

    /* could be corrupted by D_MV[*], so recompute: */
    rd_gse(&dom, D_MP0, rags, rage);
    /* extend bulk borders if density is outside */
    for (i=0; i<dim; i++) {
      rd_gse(&dom,D_MV[i],raVgs,raVge);
      for (j=0; j<dim; j++) {
	if (j == i) {
	  if ( raVgs[j] < rags[j] ) rags[j] = raVgs[j];
	  if ( raVge[j] >=rage[j] ) rage[j] = raVge[j] + 1;
	}
	else {
	  if ( raVgs[j] < rags[j] ) rags[j] = raVgs[j];
	  if ( raVge[j] > rage[j] ) rage[j] = raVge[j];
	}			
      }  
    }
    err = ra_create(&Btmp, dim, rags, rage);
    if ( err ) {
      fprintf(stream,
	      "ERROR: sg_readmedia from rsfread - cannot allocate tmp array.\n");
      return err;
    }
    ra_size(&Btmp, ran);
	  
    buoykey=NULL;
    rhokey=NULL;

    if (!ps_flcstring(par,"buoyancy",&buoykey)) {
      err=rsfread(Btmp._s,rags,ran,buoykey,1,stream,panelindex);
      userfree_(buoykey);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia from rsfread - buoykey (single) = %s\n",buoykey);
	ra_destroy(&Btmp);
	return err;
      }

      /* buoyancy check */
      vmax=Btmp._s0[0];
      vmin=Btmp._s0[0];
      for (j=1;j<ntot;j++) {
	vmax=iwave_max(vmax,Btmp._s0[j]);
	vmin=iwave_min(vmin,Btmp._s0[j]);
      }
      err = (vmin < bmin) || (vmax > bmax);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia - buoyancyy field read from file violates bounds [%e, %e]\n",bmin,bmax);
	fprintf(stream,"min buoyancy = %e max buoyancy = %e\n",vmin,vmax);
	
	return E_OTHER;
      }

      /* shift along coord axes, average */
      for (i = 0; i < dim; i++) {
	rd_gse(&dom, D_MV[i], raVgs, raVge);
	for ( j = dim; j < RARR_MAX_NDIM; ++j ) 
	  raVgs[j] = raVge[j] = 0; /* extra dimensions */

#if RARR_MAX_NDIM > 2
	for ( ip[2] = raVgs[2]; ip[2] <= raVge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
	  for ( ip[1] = raVgs[1]; ip[1] <= raVge[1]; ++ip[1] )
#endif
	    for ( ip[0] = raVgs[0]; ip[0] <= raVge[0]; ++ip[0] ) {
	      q = ra_gget(&Btmp, ip);
	      ip[i] += 1;
	      q += ra_gget(&Btmp, ip);
	      ip[i] -= 1;
	      
	      rd_gset(&dom, D_MV[i], ip, q * 0.5);
	    }
      }	  
    }
    
    else if ((!ps_flcstring(par,"density",&rhokey)) || 
	     (!ps_flcstring(par,"rho",&rhokey))) {
      //      fprintf(stream,"sg_readmedia -> rsfread dens\n");
      //      fflush(stream);
      err=rsfread(Btmp._s,rags,ran,rhokey,1,stream,panelindex);
      userfree_(rhokey);
      //      fprintf(stream,"sg_readmedia <- rsfread dens\n");
      //      fflush(stream);
      if (err) {
	fprintf(stream,"ERROR: sg_readmedia from ra_rsfread - rhokey (single) = %s\n",rhokey);
	ra_destroy(&Btmp);
	return err;
      }

      /* density check */
      vmax=Btmp._s0[0];
      vmin=Btmp._s0[0];
      for (j=1;j<ntot;j++) {
	vmax=iwave_max(vmax,Btmp._s0[j]);
	vmin=iwave_min(vmin,Btmp._s0[j]);
      }
      err = (vmin < dmin) || (vmax > dmax);
      if (err) {
	fprintf(stream,
		"ERROR: sg_readmedia - density field read from file violates bounds [%e, %e]\n",dmin,dmax);
	fprintf(stream,"min density = %e max density = %e\n",vmin,vmax);
	fflush(stream);
	return E_OTHER;
      }

      /* shift along coord axes, average */
      for (i = 0; i < dim; i++) {
	rd_gse(&dom, D_MV[i], raVgs, raVge);
	for ( j = dim; j < RARR_MAX_NDIM; ++j ) 
	  raVgs[j] = raVge[j] = 0; /* extra dimensions */

#if RARR_MAX_NDIM > 2
	for ( ip[2] = raVgs[2]; ip[2] <= raVge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
	  for ( ip[1] = raVgs[1]; ip[1] <= raVge[1]; ++ip[1] )
#endif
	    for ( ip[0] = raVgs[0]; ip[0] <= raVge[0]; ++ip[0] ) {
	      q = ra_gget(&Btmp, ip);
	      ip[i] += 1;
	      q += ra_gget(&Btmp, ip);
	      ip[i] -= 1;
	      
	      rd_gset(&dom, D_MV[i], ip, 2.0 / q);
	    }
      }	  
    }
    else {
      /* last resort: default value */
      for (i=0;i<dim;i++) {
	rd_size(&dom,D_MV[i],ran);
	m=1;
	for (j=0;j<dim;j++) m*=ran[j];
	for (j=0;j<m;j++) dom._s[D_MV[i]]._s0[j]=refbuoy;
      }
    }
  

    /* create bulk mod if velo provided 
       
    if buoykey, then buoyancy has been read into Btmp
    
    if rhokey,  then density  has been read into Btmp
    
    if neither, then nothing  has been read into Btmp and
    don't use it.
    
    */
    /* NOTE POSSIBILITY OF ZERODIVIDE - SHOULD BE PROTECTED */
    
    if (veloflag) {
      /* default density value */
      b = REAL_ONE/refbuoy;

      rd_gse(&dom, D_MP0, raPgs, raPge);
      for ( j = dim; j < RARR_MAX_NDIM; ++j ) 
	raPgs[j] = raPge[j] = 0; /* extra dimensions */
#if RARR_MAX_NDIM > 2
      for ( ip[2] = raPgs[2]; ip[2] <= raPge[2]; ++ip[2] ) 
#endif
#if RARR_MAX_NDIM > 1
	for ( ip[1] = raPgs[1]; ip[1] <= raPge[1]; ++ip[1] ) 
#endif
	  for ( ip[0] = raPgs[0]; ip[0] <= raPge[0]; ++ip[0] ) {
	    q = rd_gget(&dom, D_MP0, ip);
	    if (buoykey || rhokey) {
	      b = ra_gget(&Btmp, ip);
	      if (buoykey) b=REAL_ONE/b;
	    }
	    rd_gset(&dom, D_MP0, ip, q * q * b);
	  }
    }
  			  
    ra_destroy(&Btmp);
  
  }

  return err;
}

/*----------------------------------------------------------------------------*/
/*   
HERE I POPULATE ETA ARRAYS FROM THE MASTER INPUT FILE: 
*/
static int sgn_setetas(PARARRAY *pars, FILE *stream, IMODEL *model) {

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
  if ( ps_flreal(*pars, "cmax", &cmax) ) 
    fprintf(stream,"NOTE: sgn_setetas - using default cmax = %e\n",cmax);

  /* read pml amplitude ----------------------------------------------------*/
  pmlampl = 1.5 * log(1000.0);
  if ( ps_flreal(*pars,"npml_ampl", &pmlampl) )
    fprintf(stream, "NOTE. Cannot read npml_ampl. Default = %g.\n",  pmlampl);
  else
    fprintf(stream, "NOTE. Eta multiplier nplm_ampl = %g.\n", pmlampl);
    
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
	      s = eta0 * fabs(s * s * s);
            }
	  else if ( s > v1 )
            {
	      s = (s - v1) / L1;
	      s = eta1 * fabs(s * s * s);
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
	      s = eta0 * fabs(s * s * s);
            }
	  else if ( s > v1 )
            {
	      s = (s - v1) / L1;
	      s = eta1 * fabs(s * s * s);
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
int asg_readmedia(PARARRAY * pars, 
		  FILE * stream,
		  IMODEL * model,
		  int panelindex) {

  int err=0;

  /* read main medium files */
  err=sg_readmedia(model->ld_a,
		   stream,
		   *pars,
		   panelindex);

  if (!err) err=sgn_setetas(pars,stream,model);
  return err;
}

/* helper function - not in header file */
int asg_readpmlgrid(IPNT nl, IPNT nr, 
		    RPNT dx, size_t ndim, 
		    PARARRAY par, FILE * stream) {
  int idim;
  ireal fpeak = FPEAK_DEF;
  ireal cmax  = CMAX_DEF;
  ireal tmp;
  char key[30];
  size_t kl=4;

  if (ps_flreal(par,"fpeak",&fpeak))
    fprintf(stream,"NOTE: sgn_readpmlgrid - using default fpeak = %e\n",fpeak);

  if (ps_flreal(par,"cmax",&cmax)) 
    fprintf(stream,"NOTE: sgn_readpmlgrid - using default cmax = %e\n",cmax);

  for (idim=0; idim<ndim; idim++ ) {
	  if (dx[idim]<REAL_EPS) {
		  fprintf(stream, "dx[%d] = %f\n", idim, dx[idim]);
		  return E_BADINPUT;
	  }
    snprintf(key,kl,"nl%d",idim+1);
    tmp=0.0;
    ps_flreal(par, key, &tmp);

    nl[idim]=iwave_max(0,(int)((ceil)(tmp*cmax/(fpeak*dx[idim]))));
    snprintf(key,kl,"nr%d",idim+1);
    tmp=0.0;
    ps_flreal(par,key,&tmp);
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
/* old signature
int sgn_readgrid(grid * g,
		 IPNT nl,
		 IPNT nr,
		 FILE *stream, 
		 PARARRAY *pars) {
*/
int asg_readgrid(PARARRAY * pars, 
		 FILE *stream, 
		 IMODEL * model) {

  // workspace
  int err=0;
  RPNT d;
	
  IASN(model->nls, IPNT_0);  /* default left PML width */
  IASN(model->nrs, IPNT_0);  /* default right PML width */
  
  /* read physical grid parameters from bulk modulus or velocity
     grid, or from parameter file */
  if (err=init_acoustic_geom_par(&(model->g),*pars,stream)) {
    fprintf(stream,"Error: could not read grid parameters\n");
    return err;
  }

  /* extract grid steps */
  get_d(d,model->g);

  /* extract pml sizes from parameter files */
  err=asg_readpmlgrid(model->nls,model->nrs,d, (model->g).dim, *pars,stream);
  if (err) {
    fprintf(stream,"Error: could not read pml grid parameters\n");
    return err;
  }
  
  return err;
}

