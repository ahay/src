#include "utils.h"
#include "iwave.h"
#include "parser.h"
//#include "model_action.h"
#ifdef _OPENMP
#include "omp.h"
#endif


//#define VERBOSE

/* drivers should self-doc:

   #include "<modelname>_selfdoc.h"   
   xargc=argc;
   xargv=argv;
   requestdoc(1);
*/

/** IWAVE constructor
 */
// initializaion function now initializes FD_MODEL member in IMODEL
// functionally equivalent to template param or policy
// WWS 23.10.10
int iwave_construct(IWAVE * state, 
		    PARARRAY * pars,
		    FILE * stream,
		    int (*minit)(PARARRAY *, 
				 FILE *,
				 IMODEL *)) {
#ifdef VERBOSE
  fprintf(stream,"NOTE: in iwave constructor\n");
  ps_printall(*pars,stream);
  fflush(stream);
#endif
  int err=0;
  FD_MODEL * fdm = NULL;

  /* PRECONDITIONS: 
   * basic parallel environment set up, including all comms
   * output stream opened, parameter array initialized
   * minit valid function pointer for model default constructor.
   */

   /********************** BEGIN PINFO **********************************/

  /* initialize parallel info structure - part that depends only on 
     MPI_Cart_comm, not on physical grid
   */
#ifdef VERBOSE
  fprintf(stream,"initpinfo\n");
  fflush(stream);
#endif
  err=initpinfo(&(state->pinfo),stream);
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from initpinfo, err=%d\n",err);
    return err;
  }

   /********************** END PINFO **********************************/

  /* default construction of model component of state 
   */

#ifdef VERBOSE
  fprintf(stream,"im_construct\n");
  fflush(stream);
#endif
  err=im_construct(&(state->model));
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from im_construct\n");
    return err;
  }

  /* FD_MODEL initialization - should (1) dynamically allocate an
   * FD_MODEL object, (2) set pointers to FD_MODEL member functions as
   * listed in FD_MODEL, and (3) initialize void * specs data member
   * of IMODEL argument to point to FD_MODEL object. minit is should
   * conform to the fd_model_init interface in FD_MODEL. In a C++
   * implementation, a class or struct containing this and the other
   * functions could be passed as a template parameter.
   */
#ifdef VERBOSE
  fprintf(stream,"minit\n");
  fflush(stream);
#endif
  err=minit(pars,stream,&(state->model));
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from fd_model_init\n");
    return err;
  }

  /* post-construction initialization of (state->model):
   * - construct primal physical grid
   * - complete construction of cartesian grid
   * - initialization of simulator - uses IMODEL interface mcrea,
   *   which should amount to post-construction initialization..
   */

  /* extract FD_MODEL member - stored as void * in IMODEL */
  fdm = (FD_MODEL *)((state->model).specs);

  /* initialize primal grid and grid padding params */
  if ( (err=fdm->readgrid(pars, stream, &(state->model))) ) {
    fprintf(stream,"ERROR: iwave_construct from readgrid\n");
    return err;
  }

#ifdef VERBOSE
  fprintf(stream,"readparpars\n");
  fflush(stream);
#endif

  /* read misc. flags */
  readparpars(pars,
	      &(state->stats),
	      &(state->nopts),
	      &(state->printact),
	      stream);

  /* initialize data exchange info depening on physical grid dimn */

#ifdef VERBOSE
  fprintf(stream,"initexch \n");
  fflush(stream);
#endif

  err=initexch(&(state->pinfo),(state->model).g.dim,stream);
  if (err) {
    if (err==E_NOTINGRID) 
      fprintf(stream,"NOTE: Internal: create parallel returned Not-In-Grid\n");    
    else
      fprintf(stream,"ERROR: Internal create parallel error\n");

    return err;
  }

  /* post-construction initialization of (state->model) via "virtual" interface */
#ifdef VERBOSE
  fprintf(stream,"mcrea\n");
  fflush(stream);
#endif

  err=fd_modelcrea((state->pinfo).cdims,
		   (state->pinfo).crank,
		   pars,
		   stream,
		   &(state->model));
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from model create function\n");
    return err;
  }

#ifdef VERBOSE
  fprintf(stream,"setrecvexchange\n");
  fflush(stream);
#endif

  /* create recv exchange structures */
  err=setrecvexchange(&(state->model),
		      &(state->pinfo),
		      stream);
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from setrecvexchange\n");
    return err;
  }

#ifdef VERBOSE
  fprintf(stream,"prepexch\n");
  fflush(stream);
#endif

  /* initialize send/receive requests - MPI exec only */
#ifdef IWAVE_USE_MPI
  prepexch(&(state->pinfo),
	   &(state->model),
	   state->stats,
	   &(state->wt0));
#endif

#ifdef VERBOSE
  fprintf(stream,"exit iwave_construct err=%d\n",err);
  fflush(stream);
#endif

  return err;
}

/** C-style write function */
void iwave_printf(IWAVE * state, PARARRAY * pars, FILE * stream) {
  if (dump_pi(&(state->pinfo),stream)) {
    fprintf(stream,
	    "NOTE: failed to print PARALLELINFO\n");
  }
  dump_ac(pars,&(state->model),stream);
  dump_rs(pars,&(state->model),stream,0);
  dump_rs(pars,&(state->model),stream,1);
}

void iwave_destroy(IWAVE * state) {
  /* free memory allocated for workspace */
  im_destroy(&(state->model));
  /* free parallel descriptors */
  //  destroyparallel(&(state->pinfo),0); rename 10.04.11
  destroypinfo(&(state->pinfo),0); 
}

void iwave_dynamic_init(IWAVE * state,
			int istart) {
  int i;
  FD_MODEL * fdm = (FD_MODEL *)((state->model).specs);
  for (i=0;i<(state->model).ld_a.narr; i++) 
    if (fdm->isdyn(i)) 
      ra_zero(&((state->model).ld_a._s[i]));
  /* transfer start time to state */
  (state->model).tsind.it=istart;
  (state->model).tsind.iv=0;

}

int iwave_static_init(IWAVE * state, PARARRAY * par, FILE * stream, int panelindex) {
  int err=0;
  err=fd_mread(stream, 
	       &(state->model),
	       par,
	       panelindex);
  //  fprintf(stderr,"returning from iwave_static_init: ndim=%ld\n",(state->model).g.dim);
  return err;
}

int iwave_run(IWAVE * state, FILE * stream) {

  int err=0;
  int iv;
  int it;
  int ia;
  int iact;

  FD_MODEL * fdm = (FD_MODEL *)((state->model).specs);

  /* workspace needed only in MPI exec */
#ifdef IWAVE_USE_MPI
  int i;
  MPI_Status status;  
  int tmpdest, tmpsource;
  MPI_Datatype tmpdest_dt, tmpsource_dt;
  int tmpdest_val, tmpsource_val;
  void *tmpdest_buf, *tmpsource_buf;
  double time;
  time = 0.0; /* To avoid "unitialized variable" warning */
#endif

  it = (state->model).tsind.it;
  iv = (state->model).tsind.iv;

  ia = (state->model).tsinfo.pairs[iv].arr;
  iact = (state->model).tsinfo.pairs[iv].action;
    
  if ( state->printact > 1) 
  {
    fprintf(stream,"******************************************************************\n");
    fprintf(stream,"*** iwave_run: it=%d iv=%d iarr=%d iact=%d\n", it, iv, ia, iact);
    fprintf(stream,"******************************************************************\n");
  }
  
#ifdef IWAVE_USE_MPI
  if ( state->stats > 1 ) time = MPI_Wtime();
#endif
    
  switch (iact)
  {        
    case ACTION_EXCHANGE:
      if ( state->printact > 1 ) 
      {
        fprintf(stream,"  EXCHANGE array %2d", ia); 
        fflush(stream);
      }
#ifdef IWAVE_USE_MPI
      for ( i = 0; i < (state->model).nnei; ++i ) 
      {
        if ( (state->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) 
        {  
          tmpdest = (state->pinfo).sranks[i];
          tmpdest_buf = (state->pinfo).seinfo[ia][i].buf;
          tmpdest_dt = (state->pinfo).seinfo[ia][i].type;
        }
        else 
        {
          tmpdest = MPI_PROC_NULL;
          tmpdest_buf = &tmpdest_val;
          tmpdest_dt = MPI_INT;
        } 
        if ( (state->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) 
        {
          tmpsource = (state->pinfo).rranks[i];
          tmpsource_buf = (state->pinfo).reinfo[ia][i].buf;
          tmpsource_dt = (state->pinfo).reinfo[ia][i].type;
        } 
        else 
        {
          tmpsource = MPI_PROC_NULL;
          tmpsource_buf = &tmpsource_val;
          tmpsource_dt = MPI_INT;
        }
        /*fprintf(stream, "    i = %d, sending to wrk = %d, receiving from wrk = %d, [NULL = %d]\n", 
          i, tmpdest, tmpsource, MPI_PROC_NULL); fflush(stream);*/
          err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
                              tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
                              (state->pinfo).ccomm, &status);
        if ( err != MPI_SUCCESS ) 
        {
          fprintf(stream, 
            "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, it=%d, arr=%d. ABORT.\n", 
            err, i, it, ia);
          return err;
        }
      }
#else 
      if ( state->printact > 1 ) 
	{
	  fprintf(stream," --- no-op");
	}
#endif
      if ( state->printact > 1 )  
      { 
        fprintf(stream," - OK    \n "); 
        fflush(stream); 
      }
      break;
        
    case ACTION_COMPUTE:
      
      if ( state->printact > 5 ) {
	fprintf(stream,"\n------ iwave:: run Field Update iarr = %d it = %d -------------\n",ia,it);
	fprintf(stream,"iwave::run it=%d arr=%d BEFORE STEP ld_a\n",it,ia);
	rd_print(&((state->model).ld_a), ia, stream);
      }


      if ( state->printact > 1 ) 
      {
        fprintf(stream,"  iwave_run: COMPUTE, iv = %2d", iv); 
        fflush(stream);
      }
      err =  fdm->tsf
	(&((state->model).ld_c), 
	 ia, 
	 fdm->fdpars);
      if ( err > 0 ) 
	{
	  fprintf(stream, 
		  "ERROR. Internal: timestep error #%d at it=%d, arr=%d. ABORT.\n", 
		  err, it, ia);
	  return err;
	}

      if ( !state->nopts ) {
	//	if ( state->printact > 1 ) 
	//	  {
	//	    fprintf(stream,"iwave::run - COMPUTE, arr = %2d, post-timestep actions for prev step \n", ia); 
	//	    fflush(stream);
	//	  }
	// 09.11.10: WWS, after XW: pars for pts always = model
	// so - this is a poorly designed function!!!!
	// fixed 13.02.11 WWS
	// major changes of 03.12: this step deprecated - will be removed in next release
      	err = fdm->ptsf(ia,&(state->model));
      	if ( err > 0 ) {
      	  fprintf(stream, "ERROR iwave::run. Internal: post timestep error #%d at it=%d, arr=%d. ABORT.\n", err, it, ia);
	  return err;
      	}
      }

      if (state->printact > 5 ) {
	fprintf(stream,"\n iwave::run it=%d arr=%d AFTER STEP ld_a\n",it,ia);
	rd_print(&((state->model).ld_a),ia,stream);
      }
      if ( state->printact > 1 ) 
      { 
        fprintf(stream," - OK \n    "); 
	fprintf(stream,"*** giwave_dmod: exit action switch, it=%d iv=%d iarr=%d iact=%d++++++++ \n",it,iv,ia,iact);
        fflush(stream); 
      }
      break;
                
  default:
    fprintf(stream,
	    "ERROR. Internal: unknown action %d. ABORT\n", iact);
    return E_BADINPUT;
    }
  
#ifdef IWAVE_USE_MPI
  if ( state->stats > 1 ) 
  {
    time = MPI_Wtime() - time;
    if ((state->stat_times)[ia][iact][0] > time) 
    {
      (state->stat_times)[ia][iact][0] = time;
      (state->stat_calls)[ia][iact][0] = (state->stat_calls)[ia][iact][2];
    }
    if ((state->stat_times)[ia][iact][1] < time) 
    {
      (state->stat_times)[ia][iact][1] = time;
      (state->stat_calls)[ia][iact][1] = (state->stat_calls)[ia][iact][2];
    }
    (state->stat_times)[ia][iact][2] += time;
    (state->stat_calls)[ia][iact][2]++;
  }
  if ( state->printact > 1 ) {
    if ( state->stats > 1 ) fprintf(stream, "   time = %9.2e\n", time);
    else fprintf(stream, "\n");
  }
#endif
    
  /* move to driver - permits inversion code to override 
    next_step(&((state->model).tsind),(state->model).tsinfo);
  */
  if ( state->printact > 0 ) fflush(stream);
  
  return err;
}

