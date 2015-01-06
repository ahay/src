#include "utils.h"
#include "iwave.h"
#include "parser.h"
#ifdef _OPENMP
#include "omp.h"
#endif

/* #define IWAVE_VERBOSE */

/* drivers should self-doc:

   #include "<modelname>_selfdoc.h"   
   xargc=argc;
   xargv=argv;
   requestdoc(1);
*/

/** IWAVE constructor
 */
/* initializaion function now initializes FD_MODEL member in IMODEL */
/* functionally equivalent to template param or policy */
/* WWS 23.10.10 */
int iwave_construct(IWAVE * state, 
		    PARARRAY * pars,
		    FILE * stream,
		    IWaveInfo const & ic) {

#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: in iwave constructor\n");
  ps_printall(*pars,stream);
  fflush(stream);
#endif
  int err=0;

  /* PRECONDITIONS: 
   * basic parallel environment set up, including all comms
   * output stream opened, parameter array initialized
   * minit valid function pointer for model default constructor.
   */

  /********************** BEGIN PINFO **********************************/

  /* initialize parallel info structure - part that depends only on 
     MPI_Cart_comm, not on physical grid
  */
#ifdef IWAVE_VERBOSE
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

#ifdef IWAVE_VERBOSE
  fprintf(stream,"im_construct\n");
  fflush(stream);
#endif
  err=im_construct(&(state->model));
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from im_construct\n");
    return err;
  }

  /* post-construction initialization of (state->model):
   * - construct primal physical grid
   * - complete construction of cartesian grid
   * - initialization of simulator - uses IMODEL interface mcrea,
   *   which should amount to post-construction initialization..
   */

  /* initialize primal grid and grid padding params */
  if ( (err=fd_readgrid(pars, stream, &(state->model), ic.get_iwave_iokeys()[0].keyword)) ) {
    fprintf(stream,"ERROR: iwave_construct from readgrid\n");
    return err;
  }

  //  err = fdm->readtimegrid(pars, stream, &(state->model));
  //  err = fdm->readtimegrid(pars, stream, (state->model).g, ((state->model).tsind).dt);
  err = ic.get_timegrid()(pars, stream, (state->model).g, ((state->model).tsind).dt);
  if (err) {
    fprintf(stream, "Error: iwave from get_timegrid, err=%d\n",err);
    return err;
  }
  if (((state->model).tsind).dt <= REAL_ZERO) {
    fprintf(stream, "Error: iwave_construct\n");
    fprintf(stream, "  bad input: wrong time step dt=%g\n", 
	    ((state->model).tsind).dt);
    return E_BADINPUT;
  }
#ifdef VERBOSE
  fprintf(stderr, "dt = %g\n", ((state->model).tsind).dt);
#endif
  /* FD_MODEL initialization - should (1) dynamically allocate an
   * FD_MODEL object, (2) set pointers to FD_MODEL member functions as
   * listed in FD_MODEL, and (3) initialize void * specs data member
   * of IMODEL argument to point to FD_MODEL object. minit is should
   * conform to the fd_model_init interface in FD_MODEL. In a C++
   * implementation, a class or struct containing this and the other
   * functions could be passed as a template parameter.
   */
#ifdef IWAVE_VERBOSE
  fprintf(stream,"minit\n");
  fflush(stream);
#endif

  err=ic.get_minit()(*pars,stream,state->model);
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from fd_model_init\n");
    return err;
  }

  // cerr<<"numsubsteps\n";
  /* transfer number of substeps to TSINDEX member of IMODEL */
  (state->model).tsind.niv=fd_numsubsteps(ic);

#ifdef IWAVE_VERBOSE
  fprintf(stream,"readparpars\n");
  fflush(stream);
#endif

  // cerr<<"readparpars\n";

  /* read misc. flags */
  readparpars(pars,
	      &(state->stats),
	      &(state->nopts),
	      &(state->printact),
	      stream);

  /* initialize data exchange info depening on physical grid dimn */

#ifdef IWAVE_VERBOSE
  fprintf(stream,"initexch \n");
  fflush(stream);
#endif

  // cerr<<"initexch\n";
  err=initexch(&(state->pinfo),(state->model).g.dim,stream);
  if (err) {
    if (err==E_NOTINGRID) {
#ifdef IWAVE_VERBOSE
      fprintf(stream,"NOTE: Internal: create parallel returned Not-In-Grid\n");    
#endif
    }
    else
      fprintf(stream,"ERROR: Internal create parallel error\n");
      
    return err;
  }

  /* post-construction initialization of (state->model) via "virtual" interface */
#ifdef IWAVE_VERBOSE
  fprintf(stream,"mcrea\n");
  fflush(stream);
#endif

  
  // cerr<<"before modelcrea\n";
  err=fd_modelcrea((state->pinfo).cdims,
		   (state->pinfo).crank,
		   pars,
		   stream,
		   &(state->model),
		   ic);
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from model create function\n");
    return err;
  }

  // cerr<<"return from modelcrea"<<endl;
#ifdef IWAVE_VERBOSE
  fprintf(stream,"setrecvexchange\n");
  fflush(stream);
#endif

  /* create recv exchange structures */
  err=setrecvexchange(&(state->model),
		      &(state->pinfo),
		      stream,ic);
  if (err) {
    fprintf(stream,"ERROR: iwave_construct from setrecvexchange\n");
    return err;
  }

#ifdef IWAVE_VERBOSE
  fprintf(stream,"prepexch\n");
  fflush(stream);
#endif

  /* initialize send/receive requests - MPI exec only */
#ifdef IWAVE_USE_MPI
  prepexch(&(state->pinfo),
	   &(state->model),
	   ic); 
#endif

  /* dump domain info if requested */
  iwave_printf(state,pars,stream);

#ifdef IWAVE_VERBOSE
  fprintf(stream,"exit iwave_construct err=%d\n",err);
  fflush(stream);
#endif

  return err;
}

/** C-style write function */
void iwave_printf(IWAVE * state, PARARRAY * pars, FILE * stream) {
  dump_pi(pars,&(state->pinfo),stream);
  dump_ac(pars,&(state->model),stream);
  dump_rs(pars,&(state->model),stream,0);
  dump_rs(pars,&(state->model),stream,1);
}

void iwave_destroy(IWAVE * state, FD_MODELDEST mdest) {
  /* free memory allocated for workspace */
  im_destroy(&(state->model),mdest);
  /* free parallel descriptors */
  destroypinfo(&(state->pinfo),0); 
}

/* replaced ra_zero with ra_a_zero - WWS 05.09.12 */

void iwave_dynamic_init(IWAVE * state,
			int istart, IWaveInfo const & ic) {
  int i;
  for (i=0;i<(state->model).ld_a.narr; i++) 
    //	if (isdyn(fdm,i)) 
    if (fd_isdyn(i,ic)) 
      ra_a_zero(&((state->model).ld_a._s[i]));
  /* transfer start time to state */
  (state->model).tsind.it=istart;
  (state->model).tsind.iv=0;

}




