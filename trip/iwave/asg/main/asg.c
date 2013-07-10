/* Driver for coustic staggered grid modeling in 1D, 2D, or
3D. Outputs pressure traces at specified sample rates, arbitrary
receiver positions, and/or movie frames of pressure and/or
velocity. Point source at arbitrary location, or array sources with
arbitrarily many traces. Uses staggered grid finite difference scheme
of order 2 in time and 2k in space, k=1, 2, or 4, derived from
pressure-velocity form of acoustodynamics. Either pressure-free
(reflecting) or absorbing (PML) boundary conditions on boundary faces
of simulation hypercube. Optionally parallelized at both loop (domain
decomposition) and task (shot) levels.

Authors: Igor Terentyev, William W. Symes, Tetyana Vdovina, Xin Wang

Typical parameter list. May be copied, edited, and used for input: place
in file "parfile", invoke modeling command by 

sfasg par=parfile [Madagascar install] 

or 

asg.x par=parfile [standalone install]

Given values are defaults, non-optional values in corner brackets.

-------------------------------- cut here ------------------------------

INPUT DATA FOR iwave

------------------------------------------------------------------------
FD:

         order = 2          spatial half-order
           cfl = 0.75       proportion of max dt/dx
          cmin = 1.0        min permitted velocity (m/ms) - sanity check
          cmax = 4.5        max permitted velocity (m/ms) - used in dt comp
          dmin = 0.5        min permitted density (g/cm^3) - sanity check
          dmax = 3.0        max permitted density (g/cm^3) - sanity check
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

       srctype = point       source type - point or array
        source = <path>      path to input source data file (traces), 
	                     SU format - required for array source, 
			     optional for point source
       sampord = 1           spatial interp order - 0 (nr nbr) or 
                             1 (multilin)

point source options:

                             if "source" keyword not used, Ricker source 
                             generated with central freq = fpeak
         phase = zerophase   wavelet phase (Ricker option) - possible values
                             are zerophase, minimumphase, maximumphase
       refdist = 1000.0      reference distance (m) - if > 0, scale source
                             to achieve target amplitude at this distance in
                             homogeneous 3D medium.
        refamp = 1.0         target amplitude at reference (GPa) 

------------------------------------------------------------------------
Trace info:

      hdrfile  = <path>      input prototype trace file, SU format
     datafile  = <path>      output data file, SU format - headers
                             from hdrfile, trace data from modeling

------------------------------------------------------------------------
Model info:

      velocity = <path>      path to input velocity file, RSF format
       density = <path>      path to input density file, RSF format

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
#include <sgn.h>
#include <trace_term.h>
#include <pointsrc.h>
#include <sampler.h>
#include <parser.h>
#include <asg_selfdoc.h>
#include <par.h>
#include <asg_movie.h>
      
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
  POINTSRC * ptsrc=NULL;   /* acoustic pt src    */
  SAMPLER * arrsrc=NULL;   /* source array       */
  int towed=0;             /* fixed vs. towed    */
  RPNT d;                  /* workspace for i/o  */
  RPNT o;                  /* workspace for i/o  */
  MOVIE mt;                /* movie writer       */
  char * srctype;          /* string length      */
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

    /* extract source key - default = "point" */
    if (ps_flcstring(*pars,"srctype",&srctype)) {
      srctype=(char *)usermalloc_(6*sizeof(char));
      strcpy(srctype,"point");
    }

    /* construct iwave object */
#ifdef IWAVE_VERBOSE
    fprintf(stream,"iwave_construct\n");
    fflush(stream);
#endif

    err=iwave_construct(&state,pars,stream,&asg_modelinit);
    if (err) {
      fprintf(stream,"ERROR: main from iwave_construct. ABORT\n");
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

    /* assign default sampler params - still pressure only for this version */
    IASN(tindex,D_P);
    if (!((state.model.g.dim > 0) && (state.model.g.dim < RARR_MAX_NDIM+1))) {
      fprintf(stream,"ERROR: main - model dim = %d not in range [1, %d]. ABORT\n",state.model.g.dim,RARR_MAX_NDIM);
      err=E_BADINPUT;
      iwave_destroy(&state);
      abortexit(err,pars,&stream);
    }      
    for (i=0;i<RARR_MAX_NDIM;i++) tmult[i]=REAL_ONE/((ireal)(state.model.g.dim));
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
    /* construct movie object */
    err = asg_movie_construct(&mt, stream);
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
      err = iwave_static_init(&state,pars,stream,
				trace.t.tg.xrec,trace.t.tg.first);
      if (err) {  
	fprintf(stream,"ERROR: main from iwave_static_init \n xrec = %d, err = %d. ABORT\n", 
		trace.t.tg.xrec, err);
	abortexit(err,pars,&stream);
      }

      if (err) {  
	fprintf(stream,"ERROR: main from iwave_static_init \n xrec = %d, err = %d. ABORT\n", 
		trace.t.tg.xrec, err);
	abortexit(err,pars,&stream);
      }
      
      /* sampler_init should succeed in record range [first,last) */

      /*      fprint_grid(stderr,(state.model).g); */
      if (sampler_init(&trace,&(state.model),pars,stream)) {
	fprintf(stream,"ERROR: main from sampler_init. ABORT\n");
	abortexit(err,pars,&stream);
      }

      fprintf(stderr,"IWAVE::asg rkw=%d rk=%d isrc=%d\n",
              retrieveGlobalRank(),retrieveRank(),trace.t.tg.irec);

      if (dump_term) sampler_fprint(&trace, stream);
    
      /* pointsrc_init: initializes source wavelet with proper
       * calibration, determines start time of simulation 
       */
#ifdef IWAVE_VERBOSE
      fprintf(stream,"initialize source\n");
      fflush(stream);
#endif	
      if (!strcmp(srctype,"point")) {
	if (!(ptsrc=(POINTSRC *)usermalloc_(sizeof(POINTSRC)))) {
	  fprintf(stream,"ERROR: failed to allocate point source. ABORT\n");
	  abortexit(err,pars,&stream);
	}
	err=pointsrc_init(ptsrc,&(state.model),pars,&(trace.t.tg),stream);
	if (err) {
	  fprintf(stream,"ERROR: main from pointsrc_init. ABORT\n");
	  abortexit(err,pars,&stream);
	}
	if (dump_term) pointsrc_fprint(ptsrc,stream);
	istart=ptsrc->istart;
      }
      else if (!strcmp(srctype,"array")) {
	if (!(arrsrc=(SAMPLER *)usermalloc_(sizeof(SAMPLER)))) {
	  fprintf(stream,"ERROR: failed to allocate array source. ABORT\n");
	  abortexit(err,pars,&stream);
	}
	ps_flint(*pars,"towed",&towed);
	if (!towed) {
	  RASN(scoord,RPNT_0);
	}
	
	IASN(sindex,D_P);
	for (i=0;i<RARR_MAX_NDIM;i++) smult[i]=REAL_ONE;

	err=sampler_construct(arrsrc,pars,sindex,smult,scoord,-1,srckey,srckey,stream);
	if (err) {
	  fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
	  abortexit(err,pars,&stream);
	}
	err=sampler_init(arrsrc,&(state.model),pars,stream);
	if (err) {
	  fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
	  abortexit(err,pars,&stream);
	}
	istart=(arrsrc->t).istart;
	if (dump_term) sampler_fprint(arrsrc,stream);
      }
      else {
	if (srctype) fprintf(stream,"ERROR: unknown source option = %s\n",srctype);
	else fprintf(stream,"ERROR: unknown source option\n");
	abortexit(err,pars,&stream);
      }
      
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
      
	if (ptsrc) {
	  err=pointsrc_run(ptsrc,&(state.model));
	  if (err) {
	    fprintf(stream,"ERROR: main from pointsrc_run. ABORT\n");
	    abortexit(err,pars,&stream);
	  }
	}
	if (arrsrc) {
	  err=sampler_run(arrsrc,&(state.model));
	  if (err) {
	    fprintf(stream,"ERROR: main from pointsrc_run. ABORT\n");
	    abortexit(err,pars,&stream);
	  }
	}

#ifdef VERBOSE_STEP
	fprintf(stream,"poststep - traces\n");
	fflush(stream);
#endif
      
	err=sampler_run(&trace,&(state.model));
	if (err) {
	  fprintf(stream,"ERROR: main from sampler_run. ABORT\n");
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
      if (ptsrc) {
	pointsrc_destroy(ptsrc);
	userfree_(ptsrc);
      }
      if (arrsrc) {
	sampler_destroy(arrsrc);
	userfree_(arrsrc);
      }
    }

    if (srctype) userfree_(srctype);

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

