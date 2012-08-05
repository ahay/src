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
	   &(state->model)); 
  // REMOVED TEMP 13.03.12 WWS
  //,
	   //	   state->stats,
	   //	   &(state->wt0));
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
    if (isdyn(fdm,i)) 
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
  if (!err && (state->printact > 5)) {
    fprintf(stream,"----------------- after iwave_static_init->fd_mread --------------------\n");
    rd_a_print(&((state->model).ld_a),stream);
  }
  return err;
}

int iwave_run(IWAVE * state, FILE * stream) {

  int err=0;
  int iv;
  int it;
  int ia;

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

  if ( state->printact > 0  && (state->model).tsind.iv == 0) 
  {
    fprintf(stream,"it=%d\n", (state->model).tsind.it);
  }
  else if ( state->printact > 1) 
  {
    fprintf(stream,"*************************************\n");
    fprintf(stream,"* iwave_run: it=%d iv=%d\n", (state->model).tsind.it, (state->model).tsind.iv);
    fprintf(stream,"*************************************\n");
  }
  
#ifdef IWAVE_USE_MPI
  if ( state->stats > 1 ) time = MPI_Wtime();
#endif

  /* This modified design (03.12, WWS) presumes that all fields start
     out synchronized. Therefore synchronization need only occur after
     a field update. Update functions must be local. Any number of
     fields can be updated in a single update function call. The
     FD_MODEL::update function flags those fields updated in the
     current substep, which must be synchronized.
  */    
  
  if ( state->printact > 1 ) {
    fprintf(stream,"\n------ iwave_run: COMPUTE it = %d iv = %d -------------\n",it,iv);
    fflush(stream); 
  }
  if ( state->printact > 5 ) {
    fprintf(stream,"\n------ iwave_run: before update\n");
    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
      if (fdm->update(ia,iv)) {
	fprintf(stream,"------ iarr = %d\n",ia);
	rd_print(&((state->model).ld_a), ia, stream);
      }
    }
    fflush(stream); 
  }

  err =  fdm->tsf
    (&((state->model).ld_c), 
     iv, 
     fdm->fdpars);

  if ( err > 0 ) {
    fprintf(stream, 
	    "ERROR. Internal: timestep error #%d at it=%d, iv=%d. ABORT.\n", 
	    err, it, iv);
    return err;
  }

  if ( state->printact > 5 ) {
    fprintf(stream,"\n------ iwave_run: after update\n");
    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
      if (fdm->update(ia,iv)) {
	fprintf(stream,"------ iarr = %d\n",ia);
	rd_print(&((state->model).ld_a), ia, stream);
      }
    }
    fflush(stream); 
  }  

#ifdef IWAVE_USE_MPI
  if ( state->printact > 1 ) {
    fprintf(stream,"\n------ iwave_run: EXCHANGE it = %d iv = %d -------------\n",it,iv);
    fflush(stream); 
  }
  
  for (ia=0;ia<RDOM_MAX_NARR;ia++) {
    if (fdm->update(ia,iv)) {
  
      for ( i = 0; i < (state->model).nnei; ++i ) {
        if ( (state->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) {  
	  tmpdest = (state->pinfo).sranks[i];
	  tmpdest_buf = (state->pinfo).seinfo[ia][i].buf;
          tmpdest_dt = (state->pinfo).seinfo[ia][i].type;
        }
        else {
          tmpdest = MPI_PROC_NULL;
          tmpdest_buf = &tmpdest_val;
          tmpdest_dt = MPI_INT;
        } 
        if ( (state->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) {
          tmpsource = (state->pinfo).rranks[i];
          tmpsource_buf = (state->pinfo).reinfo[ia][i].buf;
          tmpsource_dt = (state->pinfo).reinfo[ia][i].type;
        } 
        else {
          tmpsource = MPI_PROC_NULL;
          tmpsource_buf = &tmpsource_val;
          tmpsource_dt = MPI_INT;
        }
	if (state->printact > 1) {
	  fprintf(stream, "    i = %d, sending to wrk = %d, receiving from wrk = %d, [NULL = %d]\n", 
		  i, tmpdest, tmpsource, MPI_PROC_NULL); fflush(stream);
	}
	
	err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
			   tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
			   (state->pinfo).ccomm, &status);
        if ( err != MPI_SUCCESS ) {
          fprintf(stream, 
		  "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, it=%d, arr=%d. ABORT.\n", 
		  err, i, it, ia);
          return err;
        }
      }
    }
  }
#endif

  /* REMOVED TEMP 13.03.12
     this is probably worth reviving somehow
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
  */

  if ( state->printact > 0 ) fflush(stream);
  
  return err;
}

