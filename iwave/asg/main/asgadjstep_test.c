/* driver to test time step funcs

   asgadjstep_test.c   
   driver to test time step funcs defined in asg24_2d.c 
*/

#include <sgn.h>
#include <iwave.h>
#include <parser.h>

#include <time.h>

//   #define VERBOSE

// declare some funcitons defined in other source files as local funcs only
int asg_fts2d_24p01(RDOM *dom, void *pars);
int asg_ftsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars);

int asg_fts2d_24v0(RDOM *dom, void *pars);
int asg_ftsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars);

int asg_fts2d_24v1(RDOM *dom, void *pars);
int asg_ftsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars);

int asg_isdyn(int i);


void assign_rand(RDOM * u,int _static) {
  // assign random number to allocated domain 
  /* _static = 0, static fields only*/
  /* _static = 1, dynamic fields only */
  /* _static > 1, all fields */
  // only assign random value to computational domain
  int ndim = 2;
  int j;               // dim counter
  int i;               // grid counter
  int iarr;            // array counter
  ireal shift;         // shift away from zero, or to zero mean
  int len=0;             // length of array

  srand(time(NULL));
 
  // loop over domains (hard wire!)
  // only randomize  D_P0, D_P1, D_V0, D_V1, D_MP0, D_MV0, D_MV1
  // assign random numbers (1) between -0.5 and 0.5 with zero mean, for dynamic fields; 
  // (2) between 1.0 and 2.0 for static fields.
 
  for (iarr=0;iarr<RDOM_MAX_NARR;iarr++) {
    if (iarr <= 5 || iarr ==10) {
      len=1;
      // extract dimensions
      for (j=0;j<ndim;j++) 
	len *= u->_s[iarr]._dims[j].n0; // allocated domain
      
      shift=1.0;
      if (asg_isdyn(iarr)) {
	shift = -0.5;
	if (_static) 
	  /* randomize dynamic fields  D_P0, D_V0, D_V1, D_P1*/
	  for (i=0;i<len;i++)
	    u->_s[iarr]._s0[i]=shift + ((ireal)(rand()))/((ireal)(RAND_MAX));  
      }
      else {
	if (_static==0 || _static>1)  
	  /* randomize static fields  D_MP0, D_MV0, D_MV1*/
	  for (i=0;i<len;i++)
	    u->_s[iarr]._s0[i]=shift + ((ireal)(rand()))/((ireal)(RAND_MAX));  	
      }
    }
  }
}

void assign_const(RDOM * u, float _val, int _static) {
  // assign const to allocated domain 
  /* _static = 0, static fields only */
  /* _static = 1, dynamic fields only */
  /* _static > 1, all fields */
  int ndim = 2;
  int j;               // dim counter
  int i;               // grid counter
  int iarr;            // array counter
  int len=0;             // length of array

  // loop over domains (hard wire!)
  // only D_P0, D_P1, D_V0, D_V1, D_MP0, D_MV0, D_MV1
  for (iarr=0;iarr<RDOM_MAX_NARR;iarr++) {
    if (iarr <= 5 || iarr ==10) {
      len=1;
      // extract dimensions
      for (j=0;j<ndim;j++) 
	len *= u->_s[iarr]._dims[j].n0; // allocated domain
      
      if (asg_isdyn(iarr)) {
	if (_static) 
	  for (i=0;i<len;i++) 
	    u->_s[iarr]._s0[i]=_val;
      }
      else {
	if ( (_static == 0) || (_static > 1) )
	  for (i=0;i<len;i++) 
	    u->_s[iarr]._s0[i]=_val;
      }
    }
  }
}

float comp_rd_inner(RDOM * rd1, RDOM * rd2, int _static) {
  /* _static = 0, compute inner product for static fields */
  /* _static != 0, compute inner product for dynamic fields */
  int ndim = 2;
  int j;               // dim counter
  int i;               // grid counter
  int iarr=0;            // array counter
  int len=0;             // length of array
  float val = 0.0f;

  // loop over domains(hard wire!)
  // only D_P0, D_P1, D_V0, D_V1, D_MP0, D_MV0, D_MV1
  for (iarr=0;iarr<RDOM_MAX_NARR;iarr++) {
    if (iarr <= 5 || iarr ==10) {
      len=1;
      // extract dimensions
      for (j=0;j<ndim;j++) 
	len *= rd1->_s[iarr]._dims[j].n0; // allocated domain

      if (asg_isdyn(iarr)) {
	if (_static) 
	  for (i=0;i<len;i++)
	    val += rd1->_s[iarr]._s0[i] * rd2->_s[iarr]._s0[i];
      }
      else {
	if ( _static == 0 ) 
	  for (i=0;i<len;i++)
	    val += rd1->_s[iarr]._s0[i] * rd2->_s[iarr]._s0[i];
      }
    }
  }
  return val;
}

void write_field(RDOM * u, int _static, FILE * stream) {
  // print out fields 
  /* _static = 0, static fields */
  /* _static = 1, dynamic fields */
  /* _static > 1, all fields */
  int j;               // dim counter
  int i;               // grid counter
  int iarr;            // array counter
  int nx, ny,nx0;

  fprintf(stream, "\n\n------------------------------------------\n");
  fprintf(stream, "------print fields --------\n");
  // loop over domains (hard wire!)
  // only zero  D_P0, D_P1, D_V0, D_V1, D_MP0, D_MV0, D_MV1
  for (iarr=0;iarr<RDOM_MAX_NARR;iarr++) {
    if (iarr <= 5 || iarr ==10) {
      nx0 = u->_s[iarr]._dims[0].n0;
      nx = u->_s[iarr]._dims[0].n0;
      ny = u->_s[iarr]._dims[1].n0;
      if (asg_isdyn(iarr)) {
	if (_static) {
	  fprintf(stream, "\n------dynamic field: iarr = %d, nx = %d, ny=%d -------------\n",iarr,nx,ny);
	  for (i=0;i<nx;i++){
	    for (j=0; j<ny;j++)
	      fprintf(stream," %8.6f ",u->_s[iarr]._s0[i + j* nx0]);
	    fprintf(stream,"\n");
	  }
	}
      }
      else {
	if ( (_static == 0) || (_static > 1) ){
	  fprintf(stream, "\n------static field: iarr = %d, nx = %d, ny=%d -------------\n",iarr,nx,ny);
	  for (i=0;i<nx;i++){
	    for (j=0; j<ny;j++)
	      fprintf(stream," %8.6f ",u->_s[iarr]._s0[i + j* nx0]);
	    fprintf(stream,"\n");
	  }
	}
      }
    }
  }
}


int main(int argc, char ** argv) {

  int err=0;               /* error flag         */
  FILE * stream;           /* output stream      */
  PARARRAY pars;           /* parameter array    */
  int ts = 0;
  int rk = 0;
  IWAVE rstate, dstate, astate;      /* model states  */
  float val_inner1, val_inner2, val_nmsum;
  
  // values used to scale space and time inner products
  float sp_vol = 1.0f;
  float  _dt = 1.0f;
 
  SGN_TS_PARS * _tspars = NULL;
  FD_MODEL * _fdm = NULL; 
  //  RPNT dd;

  /* initialize parallel environment, output stream, and param table */
#ifdef IWAVE_USE_MPI
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &ts);
#endif

  initparallel_global(ts);
  rk=retrieveGlobalRank();
#ifdef VERBOSE
  fprintf(stderr,"Global MPI_Comm_size = %d\n",retrieveGlobalSize());
  fprintf(stderr,"initoutstream\n");
#endif
  err=initoutstream(&stream,retrieveGlobalRank(),retrieveGlobalSize());
  if (err) {
    fprintf(stderr,"ERROR: main from initoutstream. ABORT\n");
    abortexit(err,&pars,&stream);
  }
#ifdef VERBOSE
  fprintf(stream,"readinput\n");
  fflush(stream);
#endif
  readinput(&pars,stream,argc,argv);
  
#ifdef VERBOSE
  fprintf(stream,"paramtable:\n");
  ps_printall(pars,stream);
  fflush(stream);
#endif

#ifdef VERBOSE
  fprintf(stream,"initparallel_local \n");
  fflush(stream);
#endif
  initparallel_local(pars,stream);

  /* construct iwave object */
#ifdef VERBOSE
  fprintf(stream,"iwave_construct\n");
  fflush(stream);
#endif
  err=iwave_construct(&rstate,&pars,stream,&asg_modelinit);
  err= err || iwave_construct(&dstate,&pars,stream,&asg_modelinit);
  err= err || iwave_construct(&astate,&pars,stream,&asg_modelinit);

  if (err) {
    fprintf(stream,"ERROR: main from iwave_construct. ABORT\n");
    iwave_destroy(&rstate);
    iwave_destroy(&dstate);
    iwave_destroy(&astate);
    abortexit(err,&pars,&stream);
  }

#ifdef VERBOSE
  fprintf(stream,"rstate: iwave_printf\n");
  fflush(stream);
#endif
  iwave_printf(&rstate,&pars,stream);
#ifdef VERBOSE
  fprintf(stream,"dstate: iwave_printf\n");
  fflush(stream);
#endif
  iwave_printf(&dstate,&pars,stream);
#ifdef VERBOSE
  fprintf(stream,"astate: iwave_printf\n");
  fflush(stream);
#endif
  iwave_printf(&astate,&pars,stream);


  /****** construct working spaces **********/
#ifdef VERBOSE
  fprintf(stream,"initialize static and dynamic states\n");
  fflush(stream);
#endif
  // err = iwave_static_init(&rstate,&pars,stream,0);
  /* 
     err = err || iwave_static_init(&dstate,&pars,stream,0);
     err = err || iwave_static_init(&astate,&pars,stream,0);
     iwave_dynamic_init(&rstate,0);
     iwave_dynamic_init(&dstate,0);
     iwave_dynamic_init(&astate,0);
  */
  
  /*-- initialize rstate-dynamic with random fields --*/
  assign_rand(&rstate.model.ld_c,3);
  //assign_const(&rstate.model.ld_a,1.0f,0);

  /*-- initialize dstate-dynamic/static with zero --*/
  assign_const(&dstate.model.ld_a, 0.0f, 3);
  /*-- initialize dstate-static with random fields --*/
  assign_rand(&dstate.model.ld_c,0);

  /*-- initialize astate-static/dynamic with zero --*/
  assign_const(&astate.model.ld_a,0.0f,3);
  /*-- initialize astate-dynamic with random fields --*/
  assign_rand(&astate.model.ld_c,1);

  /*-- build tspars -- */
  _fdm = (FD_MODEL *) rstate.model.specs;
  if (_fdm) {
    _tspars = (SGN_TS_PARS *) _fdm->fdpars;
  }
  else 
    fprintf(stream,"---asgstep_test: _fdm is NULL--------------\n");

  /*--- initialize SGN_TS_PARS data members with typical values -*/
  // hard-wire!!!
  _tspars->lbc[0] = 1;
  _tspars->lbc[1] = 1;
  _tspars->rbc[0]= 1;
  _tspars->rbc[1]= 1;
  ///* 
  _tspars->dt = 1.0f;
  _tspars->lam[0]=1.0f;
  _tspars->lam[1]=1.0f;
  _tspars->k=2;
  _tspars->ndim=2;
  //*/

  /* scaling factors */
  // in this case no scaling needed
  //  get_d(dd,rstate.model.g);
  //  sp_vol = dd[0]*dd[1];
  //  _dt = _tspars->dt;
  
  if(_tspars) {
    fprintf(stream,"\n\n --------------------------------------------------\n");
    fprintf(stream,"---asgstep_test: SGN_TS_PARS info--------------\n");
    fprintf(stream,"--- dt=%f, lam[0] = %f, lam[1]=%f, k=%d, ndim=%d --\n", _tspars->dt,_tspars->lam[0],_tspars->lam[1],_tspars->k,_tspars->ndim);
  }
  else 
    fprintf(stream,"---asgstep_test: _tspars is NULL--------------\n"); 
  
  fflush(stream);

  fprintf(stream,"--- asgstep_test: scaling factors _dt = %f, sp_vol =%f --------\n",_dt, sp_vol);
  
  fprintf(stream,"\n --- asgstep_test:initial rstate.static.nmsq = %e, rstate.dynamic.nmsq = %e --------\n",sp_vol*comp_rd_inner(&(rstate.model.ld_c), &(rstate.model.ld_c),0), _dt*comp_rd_inner(&(rstate.model.ld_c), &(rstate.model.ld_c), 1));

  fprintf(stream,"\n --- asgstep_test: initial dstate.static.nmsq = %e, dstate.dynamic.nmsq = %e --------\n",sp_vol * comp_rd_inner(&(dstate.model.ld_c), &(dstate.model.ld_c),0), _dt * comp_rd_inner(&(dstate.model.ld_c), &(dstate.model.ld_c), 1));

  fprintf(stream,"\n --- asgstep_test: initial astate.static.nmsq = %e, astate.dynamic.nmsq = %e --------\n",sp_vol * comp_rd_inner(&(astate.model.ld_c), &(astate.model.ld_c),0), _dt * comp_rd_inner(&(astate.model.ld_c), &(astate.model.ld_c), 1));

#ifdef VERBOSE
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n --BEFORE COMPUTATION------------------\n");
  fprintf(stream,"\n\n ---print rstate -------------------\n");
  write_field(&(rstate.model.ld_c), 3, stream);
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n ---print dstate -------------------\n");
  write_field(&(dstate.model.ld_c), 3, stream);
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n ---print astate -------------------\n");
  write_field(&(astate.model.ld_c), 3, stream);
  fprintf(stream,"\n -------------------------------------\n");
  fflush(stream);
#endif  
  /*------------ one-step comp testing ------------ */

  /*--- forward map ---*/
  /*
  //-- D_P0 --
  err = asg_fts2d_24p01(&(rstate.model.ld_c), _tspars);
  //-- D_V0 --
  err = err || asg_fts2d_24v0(&(rstate.model.ld_c), _tspars);  
  //-- D_V1 --
  err = err || asg_fts2d_24v1(&(rstate.model.ld_c), _tspars);
  // note: now rstate holds ref static and fwd dynamic fields 
  fprintf(stream,"\n -------------------------------------\n");
  fprintf(stream,"\n--- asgstep_test:final rstate.static.nmsq = %e, rstate.dynamic.nmsq = %e --------\n",sp_vol * comp_rd_inner(&(rstate.model.ld_c), &(rstate.model.ld_c),0), _dt * comp_rd_inner(&(rstate.model.ld_c), &(rstate.model.ld_c), 1));
  */

  /*--- linearized forward map ---*/
  /*-- dp_0, dp_1 --*/
  err = asg_ftsm2d_24p01(&(dstate.model.ld_c),&(rstate.model.ld_c), _tspars);
  /*-- p_0, p_1 --*/
  err = err || asg_fts2d_24p01(&(rstate.model.ld_c), _tspars);
  /*-- dv_0 --*/
  err = err || asg_ftsm2d_24v0(&(dstate.model.ld_c),&(rstate.model.ld_c), _tspars);  
  /*-- d_v1 --*/
  err = err || asg_ftsm2d_24v1(&(dstate.model.ld_c),&(rstate.model.ld_c), _tspars);
  /* note: now dstate holds dynamic pert fields */
  fprintf(stream,"\n--- asgstep_test:final dstate.static.nmsq = %e, dstate.dynamic.nmsq = %e --------\n",sp_vol*comp_rd_inner(&(dstate.model.ld_c), &(dstate.model.ld_c),0),_dt*comp_rd_inner(&(dstate.model.ld_c), &(dstate.model.ld_c), 1));

  /*--- compute inner products for adj test ---*/
  //val_inner1, val_nmsum
  val_inner1 = _dt * comp_rd_inner(&(dstate.model.ld_c), &(astate.model.ld_c), 1);
  val_nmsum = _dt * sqrt(comp_rd_inner(&(dstate.model.ld_c), &(dstate.model.ld_c), 1) * comp_rd_inner(&(astate.model.ld_c), &(astate.model.ld_c), 1));

  /*--- adjoint test ---*/
  /*-- D_V1 --*/
  err = err || asg_atsm2d_24v1(&(astate.model.ld_c),&(rstate.model.ld_c), _tspars);
  /*-- D_V0 --*/
  err = err || asg_atsm2d_24v0(&(astate.model.ld_c),&(rstate.model.ld_c), _tspars);  
 /*-- D_P0 --*/
  err = err || asg_atsm2d_24p01(&(astate.model.ld_c),&(rstate.model.ld_c), _tspars);
  /* note: now rstate holds static adj-pert fields */
  fprintf(stream,"\n--- asgstep_test:final astate.static.nmsq = %e, astate.dynamic.nmsq = %e --------\n",sp_vol * comp_rd_inner(&(astate.model.ld_c), &(astate.model.ld_c),0), _dt * comp_rd_inner(&(astate.model.ld_c), &(astate.model.ld_c), 1));

  /*--- compute inner products for adj test ---*/
  // val_inner2
  val_inner2 = sp_vol * comp_rd_inner(&(dstate.model.ld_c), &(astate.model.ld_c), 0);

  fprintf(stream,"\n -------------------------------------\n");
  fprintf(stream,"\n ---asgstep_test, results ------------\n");
  fprintf(stream,"--- <A x, y> = %e ------------\n",val_inner1);
  fprintf(stream,"--- <x , A^T y> = %e ------------\n",val_inner2);
  fprintf(stream,"--- |A x| |y| = %e ------------\n",val_nmsum);
  fprintf(stream,"--- ( <A x, y> - <x, A^T y> )/(|A x| |y|) = %e ------------\n",(val_inner1 - val_inner2)/val_nmsum);
  fprintf(stream,"\n -------------------------------------\n");
  fflush(stream);

#ifdef VERBOSE
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n -- AFTER ONE Born-STEP------------------\n");
  fprintf(stream,"\n\n ---print rstate -------------------\n");
  write_field(&(rstate.model.ld_c), 3, stream);
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n ---print dstate -------------------\n");
  write_field(&(dstate.model.ld_c), 3, stream);
  fprintf(stream,"\n\n -------------------------------------\n");
  fprintf(stream,"\n\n ---print astate -------------------\n");
  write_field(&(astate.model.ld_c), 3, stream);
  fprintf(stream,"\n -------------------------------------\n");
  fflush(stream);
#endif

#ifdef VERBOSE
    fprintf(stream,"iwave_destroy\n");
    fflush(stream);
#endif

    /*---------- destroy static objects and exit --------*/
    iwave_destroy(&rstate);
    iwave_destroy(&dstate);
    iwave_destroy(&astate);

#ifdef IWAVE_USE_MPI

#ifdef VERBOSE
  fprintf(stream,"MPI_Finalize\n");
  fflush(stream);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif


#ifdef VERBOSE
  fprintf(stream,"quietexit\n");
  fflush(stream);
#endif

  quietexit(&pars,&stream);

  exit(0);
}

