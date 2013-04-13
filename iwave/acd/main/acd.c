/* Driver for acoustic constant density centered difference modeling in 
2D or 3D. Outputs pressure traces at specified sample rates, arbitrary
receiver positions, and/or movie frames of pressure. Array sources with
arbitrarily many traces at arbitrary locations, receivers at arbitrary 
locations. Uses centered difference scheme
of order 2 in time and 2k in space, k=1, 2, or 4, derived from
2nd order pressure form of acoustodynamics. Either pressure-free
(reflecting) or absorbing (PML) boundary conditions on boundary faces
of simulation hypercube. Optionally parallelized at both loop (domain
decomposition) and task (shot) levels.

Authors: Muhong Zhou, William W. Symes

Typical parameter list. May be copied, edited, and used for input: place
in file "parfile", invoke modeling command by 

sfacd par=parfile [Madagascar install] 

or 

acd.x par=parfile [standalone install]

Given values are defaults, non-optional values in corner brackets.

-------------------------------- cut here ------------------------------

INPUT DATA FOR iwave

------------------------------------------------------------------------
FD:

         order = 2          spatial half-order
           cfl = 0.75       proportion of max dt/dx
          cmin = 1.0        min permitted velocity (m/ms) - sanity check
          cmax = 4.5        max permitted velocity (m/ms) - used in dt comp
         fpeak = 0.010      nominal central frequency (KHz)

------------------------------------------------------------------------
PML info: layer thickness in WL at fpeak, cmax. Set = 0.0 for reflecting
boundary

           nl1 = 0.0         axis 1 - low end
           nr1 = 0.0         axis 1 - high end
           nl2 = 0.0         axis 2 - low end
           nr2 = 0.0         axis 2 - high end
           nl3 = 0.0         axis 3 - low end
           nr3 = 0.0         axis 3 - high end

------------------------------------------------------------------------
Source info:

        source = <path>      path to input source data file (traces), 
	                     SU format
       sampord = 1           spatial interp order - 0 (nr nbr) or 
                             1 (multilin)

------------------------------------------------------------------------
Trace info:

      hdrfile  = <path>      input prototype trace file, SU format
     datafile  = <path>      output data file, SU format - headers
                             from hdrfile, trace data from modeling

------------------------------------------------------------------------
Movie info:
    moviestep  = 0.0         time between movie frames (<= 0.0, 
                             or do not include, for no movie)

------------------------------------------------------------------------
Model info:

      velocity = <path>      path to input velocity file, RSF format

------------------------------------------------------------------------
MPI info:

       mpi_np1 = 1           number of subdomains along axis 1
       mpi_np2 = 1           number of subdomains along axis 2
       mpi_np3 = 1           number of subdomains along axis 3
       partask = 1           number of shots to execute in parallel

------------------------------------------------------------------------
Output info:

      printact = 0           output at every time step
                             0 - none
                             1 - time step index
                             2 - diagnostic messages from main routines
                             > 2 - much more, useful only for serious 
                                debugging
       dump_pi = 0           dump parallel/dom. decomp info
      dump_lda = 0           dump grid data for allocated arrays
      dump_ldc = 0           dump grid data for computational arrays
      dump_ldr = 0           dump grid data for receive buffers
      dump_lds = 0           dump grid data for send buffers
     dump_term = 0           dump trace header data

-------------------------------- cut here ------------------------------
*/
/*============================================================================
 *                            BEGIN INCLUDES 
 *============================================================================*/

#include <iwave.h>
#include <acd.h>
#include <trace_term.h>
#include <sampler.h>
#include <parser.h>
#include <acd_selfdoc.h>
#include <par.h>
      
#ifdef _OPENMP
#include <omp.h>
#endif
      
#define NSTR 128
      
/* uncomment to write to the rk-dep output stream at every time step 
#define VERBOSE_STEP
 */


/*============================================================================
 *                             END INCLUDES 
 *============================================================================*/
      
/*  definitions of global variables */
						     //int xargc; char **xargv;

int main(int argc, char ** argv) {

  int err=0;               /* error flag         */
  int runflag=1;           /* run time loop if 1 */
  FILE * stream;           /* output stream      */
  PARARRAY * pars=NULL;    /* parameter array    */
  IWAVE state;             /* model state        */
  SAMPLER trace;           /* trace sampler      */
  SAMPLER * arrsrc=NULL;   /* source array       */
  int towed=0;             /* fixed vs. towed    */
  RPNT d;                  /* workspace for i/o  */
  RPNT o;                  /* workspace for i/o  */
  MOVIE mt;                /* movie writer       */
  char hdrkey[NSTR];       /* key = hdrfile      */
  char datakey[NSTR];      /* key = datafile     */
  char srckey[NSTR];       /* key = source       */
  IPNT tindex;             /* sample array index */
  RPNT tmult;              /* multiplier array   */
  IPNT sindex;             /* sample array index */
  RPNT smult;              /* multiplier array   */
  RPNT scoord;             /* source cell loc    */
  int dump_term=0;         /* trace info dump    */
  int istart=0;            /* start index        */
  int ts;                  /* thread support lvl */
  int rk;                  /* process rank       */

  int i;

  /*******************************************************************
   ****************** INITIALIZE PARALLEL ENVIRONMENT ****************
   * initialize parallel environment, output stream, and param table *
   *******************************************************************/

  xargc=argc;
  xargv=argv;

  ts=0;
#ifdef IWAVE_USE_MPI
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &ts);
#endif

  /* initparallel(ts); */
  initparallel_global(ts);

  /* self-doc if no args */
  rk=retrieveGlobalRank();
  if (rk==0) requestdoc(1); 

#ifdef IWAVE_VERBOSE
#ifdef IWAVE_USE_MPI
  if (rk==0) {
    fprintf(stderr,"Global MPI_Comm_size = %d\n",retrieveGlobalSize());
  }
#endif
  fprintf(stderr,"initoutstream\n");
#endif

  err=initoutstream(&stream,retrieveGlobalRank(),retrieveGlobalSize());
  if (err) {
    fprintf(stderr,"ERROR: main from initoutstream. ABORT\n");
    abortexit(err,pars,&stream);
  }

#ifdef IWAVE_VERBOSE
  fprintf(stream,"readinput\n");
  fflush(stream);
#endif

  err=readinput(&pars,stream,argc,argv);
  if (err) {
    fprintf(stderr,"ERROR: main from readinput. ABORT\n");
    abortexit(err,pars,&stream);
  }

#ifdef IWAVE_VERBOSE
  fprintf(stream,"paramtable:\n");
  ps_printall(*pars,stream);
  fflush(stream);
#endif

#ifdef IWAVE_VERBOSE
  fprintf(stream,"initparallel_local \n");
  fflush(stream);
#endif
  initparallel_local(*pars,stream);

  /*******************************************************************
   ****************** INITIALIZE PARALLEL ENVIRONMENT ****************
   *                         finished                                *
   *******************************************************************/

  /* check comm status - skip the rest if no comm */
#ifdef IWAVE_USE_MPI
  if (retrieveGroupID() == MPI_UNDEFINED) {
    fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
  }
  else {
#endif

    //    xargc=argc;
    /* assign key strings */
    strcpy(hdrkey,"hdrfile");
    strcpy(datakey,"datafile");
    strcpy(srckey,"source");

    /* extract trace sampler verbosity flag */
    ps_flint(*pars,"dump_term",&dump_term);

    /* construct iwave object */
#ifdef IWAVE_VERBOSE
    fprintf(stream,"iwave_construct\n");
    fflush(stream);
#endif

    err=iwave_construct(&state,pars,stream,&acd_modelinit);
    if (err) {
      fprintf(stream,"ERROR: main from iwave_construct. ABORT\n");
      iwave_destroy(&state);
      abortexit(err,pars,&stream);
    }

    /* sanity-test dimension */
    if (!((state.model.g.dim > 0) && (state.model.g.dim < RARR_MAX_NDIM+1))) {
      fprintf(stream,"ERROR: main - model dim = %d not in range [1, %d]. ABORT\n",state.model.g.dim,RARR_MAX_NDIM);
      err=E_BADINPUT;
      iwave_destroy(&state);
      abortexit(err,pars,&stream);
    }      

#ifdef IWAVE_VERBOSE
    fprintf(stream,"iwave_printf\n");
    fflush(stream);
#endif

    iwave_printf(&state,pars,stream);

#ifdef IWAVE_VERBOSE
    fprintf(stream,"data sampler\n");
    fflush(stream);
#endif

    /* assign default sampler params */
    IASN(tindex,IPNT_0);
    for (i=0;i<RARR_MAX_NDIM;i++) tmult[i]=REAL_ONE;
    RASN(scoord,RPNT_0);
  
    /* construct traceterm object */
    err=sampler_construct(&trace,pars,tindex,tmult,scoord,0,hdrkey,datakey,stream);
    if (err) {
      fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
      abortexit(err,pars,&stream);
    }
    /* initial reassignment of source coordinates */
    RASN(scoord,trace.t.tg.src[trace.t.tg.xrec]);

#ifdef IWAVE_VERBOSE
    fprintf(stream,"movie sampler\n");
    fflush(stream);
#endif
    /* add "movie1 = p" to par list - only p can be output! */
    ps_slcstring(*pars,"movie1","p");
    /* construct movie object */
    err = acd_movie_construct(&mt, stream);
    if (err) {
      fprintf(stream,"ERROR: main from movie_construct. ABORT\n");
      abortexit(err,pars,&stream);
    }

    /* postconditon: tracegeom trace.tg now equipped with 
     * - number of shot records
     * - shot coordinates for each record
     * - file offset for each record
     */

    /* NOTE: pointsrc object has trivial default construction */

#ifdef IWAVE_VERBOSE
    fprintf(stream,"loop\n");
    fflush(stream);
#endif

    /* simulation loop
     * roles of "while" conditions:
     * tg.xrec = index of next source to be computed
     * tg.irec = index of current source being computed
     * tg.nrec = index of final source 
     * sampler_init updates these via call to init_tracegeom
     * sampler_init precondition (before sim loop): xrec=irec=first
     * (index of first record in current group)
     * sampler_init postcondition: irec=xrec; xrec++;
     * so at bottom of step i (=0,1,2...) through sim loop,
     * irec=i, xrec=i+1.
     * - test xrec against last: for versions which permit incoherent 
     *   models, must preceed call to iwave_static_init,
     *   even though it is tested also in traceterm_init, because
     *   latter preceeds former.
     * - iwave_static_init: uses source coordinates to read in proper
     *   coefficient arrays, which may depend on source. computes
     *   time step satisfying stability conditions for scheme.
     * - sampler_init: reads trace geometry for shot, allocates and 
     *   initializes sampling-related arrays
     */

    while ( (trace.t.tg.xrec < trace.t.tg.last+1)) {

      /* iwave_static_init should succeed in record range [0,nrec) */
#ifdef IWAVE_VERBOSE
      fprintf(stream,"iwave_static_init\n");
      fflush(stream);
#endif
      err = iwave_static_init(&state,pars,stream,trace.t.tg.xrec,trace.t.tg.first);
      if (err) {  
	fprintf(stream,"ERROR: main from iwave_static_init \n xrec = %d, err = %d. ABORT\n", 
		trace.t.tg.xrec, err);
	abortexit(err,pars,&stream);
      }

      /* sampler_init should succeed in record range [first,last) */

      /*      fprint_grid(stderr,(state.model).g); */
#ifdef IWAVE_VERBOSE
      fprintf(stream,"sampler_init (trace)\n");
      fflush(stream);
#endif
      if (sampler_init(&trace,&(state.model),pars,stream)) {
	fprintf(stream,"ERROR: main from sampler_init. ABORT\n");
	abortexit(err,pars,&stream);
      }

      fprintf(stderr,"IWAVE::acd rkw=%d rk=%d isrc=%d\n",
              retrieveGlobalRank(),retrieveRank(),trace.t.tg.irec);

      if (dump_term) sampler_fprint(&trace, stream);
    
      /* pointsrc_init: initializes source wavelet with proper
       * calibration, determines start time of simulation 
       */
#ifdef IWAVE_VERBOSE
      fprintf(stream,"initialize source\n");
      fflush(stream);
#endif	

      if (!(arrsrc=(SAMPLER *)usermalloc_(sizeof(SAMPLER)))) {
	fprintf(stream,"ERROR: failed to allocate array source. ABORT\n");
	abortexit(err,pars,&stream);
      }
      /* assign relative coordinates (scoord) for array sources. 
	 for fixed source array, towed=0, scoord is RPNT_0.
	 for towed source, receiver coords in source input file interpreted 
         as relative to source coordinate from trace headers, assigned to
	 scoord
      */
      ps_flint(*pars,"towed",&towed);
      if (!towed) RASN(scoord,RPNT_0);
      else RASN(scoord,trace.t.tg.src[trace.t.tg.xrec]);
      
      IASN(sindex,IPNT_0);
      for (i=0;i<RARR_MAX_NDIM;i++) smult[i]=REAL_ONE;

      err=sampler_construct(arrsrc,pars,sindex,smult,scoord,-1,srckey,srckey,stream);
      if (err) {
	fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
	abortexit(err,pars,&stream);
      }
      err=sampler_init(arrsrc,&(state.model),pars,stream);
      if (err) {
	fprintf(stream,"ERROR: main from sampler_init. ABORT\n");
	abortexit(err,pars,&stream);
      }
      istart=(arrsrc->t).istart;
      if (dump_term) sampler_fprint(arrsrc,stream);
      
#ifdef IWAVE_VERBOSE
      fprintf(stream,"initialize movie\n");
      fflush(stream);
#endif
      err=movie_init(&mt, &(state.model), pars, &(trace.t.tg), stream);
      if (err) {
	fprintf(stream,"ERROR: main from movie_init. ABORT\n");
	abortexit(err,pars,&stream);
      }
      if (dump_term) movie_fprint(&mt, stream);

#ifdef IWAVE_VERBOSE
      fprintf(stream,"initialize state\n");
      fflush(stream);
#endif

      /* initialize dynamic fields */
      iwave_dynamic_init(&state,istart);

#ifdef IWAVE_VERBOSE
      fprintf(stream,"time loop\n"); 
      fflush(stream);
#endif

      ps_flint(*pars,"run",&runflag);
      
      /* time stepping */
      while ((state.model.tsind.it <= trace.t.istop) &&
	     (runflag)) {

#ifdef VERBOSE_STEP
	fprintf(stderr,"\n********* it = %d\n",state.model.tsind.it);
	fprintf(stream,"step - update dynamical fields\n");
	fflush(stream);
#endif
	err=iwave_run(&state,stream);
	if (err) {
	  fprintf(stream,"ERROR: main from iwave_run. ABORT\n");
	  abortexit(err,pars,&stream);
}	

#ifdef VERBOSE_STEP
	fprintf(stream,"poststep - source\n");
	fflush(stream);
#endif

	/* fix sampler index to access real current pressure array - 
	   this step should probably be part of the sampler_run 
	   code, currently is a hack using the movie field selector.
	*/
	for (i=0;i<RARR_MAX_NDIM;i++) arrsrc->sindex[i]=mt.iselect("p");
	err=sampler_run(arrsrc,&(state.model));
	if (err) {
	  fprintf(stream,"ERROR: main from sampler_run (array source). ABORT\n");
	  abortexit(err,pars,&stream);
	}

#ifdef VERBOSE_STEP
	fprintf(stream,"poststep - traces\n");
	fflush(stream);
#endif

      	for (i=0;i<RARR_MAX_NDIM;i++) trace.sindex[i]=mt.iselect("p");
	err=sampler_run(&trace,&(state.model));
	if (err) {
	  fprintf(stream,"ERROR: main from sampler_run (trace output). ABORT\n");
	  abortexit(err,pars,&stream);
	}

#ifdef VERBOSE_STEP
	fprintf(stream,"poststep - movie\n");
	fflush(stream);
#endif

	/* note: only zero panel recorded, even for extended model; */
	err=movie_run(&mt,&(state.model),stream,0);
	if (err) {
	  fprintf(stream,"ERROR: main from movie_run. ABORT\n");
	  abortexit(err,pars,&stream);
	}

#ifdef VERBOSE_STEP
	fprintf(stream,"update time step info\n");
	fflush(stream);
#endif
	next_step(&(state.model.tsind));
      
      }

      /*
	printtime(state.stat_times, state.stat_calls, &state, stream);
      */
      /* after time step loop, write data - this block of code should be
	 part of the operator() method of the FO built by the output
	 policy class in TSOpt.
      */
    
      get_d(d, state.model.g);
      get_o(o, state.model.g);

#ifdef IWAVE_VERBOSE
      fprintf(stream,"writetraces\n");
      fflush(stream);
#endif

      if (runflag) err=writetraces(&(trace.t.tg),d,o,stream);

      if (err) {
	fprintf(stream,"ERROR: main from tracegeom::writetraces. ABORT\n");
	abortexit(err,pars,&stream);
      }

#ifdef IWAVE_VERBOSE
      fprintf(stream,"destroy source objects\n");
      fflush(stream);
#endif
      if (arrsrc) {
	sampler_destroy(arrsrc);
	userfree_(arrsrc);
      }
    }

    /* destroy static objects and exit */

#ifdef IWAVE_VERBOSE
    fprintf(stream,"movie_destroy\n");
    fflush(stream);
#endif

    movie_destroy(&mt);

#ifdef IWAVE_VERBOSE
    fprintf(stream,"sampler_destroy\n");
    fflush(stream);
#endif

    sampler_destroy(&trace);

#ifdef IWAVE_VERBOSE
    fprintf(stream,"iwave_destroy\n");
    fflush(stream);
#endif

    iwave_destroy(&state);

    iwave_fdestroy();

#ifdef IWAVE_USE_MPI

  } /* end nontriv comm branch */

#ifdef IWAVE_VERBOSE
  fprintf(stream,"MPI_Finalize\n");
  fflush(stream);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

#ifdef IWAVE_VERBOSE
  fprintf(stream,"quietexit\n");
  fflush(stream);
#endif

  quietexit(pars,&stream);

  exit(0);
}

