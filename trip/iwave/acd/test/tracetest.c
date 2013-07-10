/* 
   tracetest.c
   William W. Symes, Oct 2012
   test driver for iwave with acoustic const density centered grid model
   after iwave.c by Igor Terentyev

   purpose: standalone test of grid memory allocation for acd package
*/
/*============================================================================
 *                            BEGIN INCLUDES 
 * ============================================================================*/

#include <iwave.h>
#include <acd.h>
#include <trace_term.h>
#include <sampler.h>
#include <parser.h>
//#include <acd_selfdoc.h>
//#include <acd_movie.h>

char * sdoc[] = {
  NULL };

#ifdef _OPENMP
#include <omp.h>
#endif

#define NSTR 128

/* uncomment to write to the rk-dep output stream at every time step 
#define VERBOSE_STEP
#define IWAVE_VERBOSE
 */

/*============================================================================
 *                             END INCLUDES 
 * ============================================================================*/

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

  // special params for test
  char * fname_pars;
  char * fname_rsf;
  char * fname_data;
  FILE * fp_pars;
  FILE * fp_rsf;
  FILE * fp_data;
  ireal * v;
  int tmp_argc;
  char ** tmp_argv;

  /*******************************************************************
   ****************** INITIALIZE PARALLEL ENVIRONMENT ****************
   * initialize parallel environment, output stream, and param table *
   *******************************************************************/

  int xargc=argc;
  char ** xargv=argv;

  ts=0;
#ifdef IWAVE_USE_MPI
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &ts);
#endif

  //initparallel(ts);
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

  // set up parfile on rk 0
  fname_pars = (char *)usermalloc_(10*sizeof(char));
  strcpy(fname_pars,"tracepar");
  fp_pars = iwave_fopen(&fname_pars,"w",NULL,stream);
  fprintf(fp_pars,"order      = 2\n");
  fprintf(fp_pars,"cfl        = 0.8\n");
  fprintf(fp_pars,"cmin       = 1.0\n");
  fprintf(fp_pars,"cmax       = 6.0\n");
  fprintf(fp_pars,"velocity   = v.rsf\n");
  fprintf(fp_pars,"hdrfile    = camhdr1.su\n");
  fprintf(fp_pars,"datafile   = data.su\n");
  fprintf(fp_pars,"mpi_np1    = 1\n");
  fprintf(fp_pars,"mpi_np2    = 1\n");
  fprintf(fp_pars,"mpi_np3    = 1\n");
  fprintf(fp_pars,"partask    = 1\n");
  fprintf(fp_pars,"dump_pi    = 1\n");
  fprintf(fp_pars,"dump_lda   = 1\n");
  fprintf(fp_pars,"dump_ldc   = 1\n");
  fprintf(fp_pars,"dump_term  = 1\n");
  fflush(fp_pars);
  iwave_fclose(fp_pars);

  fname_rsf=(char *)usermalloc_(10*sizeof(char));
  strcpy(fname_rsf,"v.rsf");
  fp_rsf = iwave_fopen(&fname_rsf,"w",NULL,stream);
  fprintf(fp_rsf,"n1=101 d1=10.0 o1=-40.0\n");
  fprintf(fp_rsf,"n2=101 d2=20.0 o2=0.0\n");
  fprintf(fp_rsf,"data_format=native_float\n");
  fprintf(fp_rsf,"data_type=velocity\n");
  fprintf(fp_rsf,"in=./v.rsf@\n");
  fflush(fp_rsf);
  iwave_fclose(fp_rsf);
  userfree_(fname_rsf);
  
  fname_data=(char *)usermalloc_(10*sizeof(char));
  strcpy(fname_data,"v.rsf@");
  fp_data = iwave_fopen(&fname_data,"w",NULL,stream);
  v = (float *)usermalloc_(101*101*sizeof(float));
  for (i=0;i<101*101;i++) {
    v[i]=2.0;
  }
  fwrite(v,sizeof(float),101*101*101,fp_data);
  fflush(fp_data);
  iwave_fclose(fp_data);
  userfree_(fname_data);

  // create header traces

  system("sunull nt=701 ntr=80 dt=0.00115 | sushw key=sx a=440 c=110 j=80 | sushw key=gx a=0 b=10 j=80 | sushw key=delrt a=0 | sushw key=selev a=-40 | sushw key=gelev a=-80 | suchw key1=offset key2=gx key3=sx a=0 b=1 c=-1 > camhdr1.su");

  tmp_argc=2;
  tmp_argv = (char **)malloc(tmp_argc*sizeof(char*));
  tmp_argv[0] = (char *)malloc((strlen(argv[0])+1)*sizeof(char));
  strcpy(tmp_argv[0],argv[0]);
  tmp_argv[1] = (char *)malloc(128*sizeof(char));
  strcpy(tmp_argv[1],"par=tracepar");

  err=readinput(&pars,stream,tmp_argc,tmp_argv);
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

    err=iwave_construct(&state,pars,stream,&acd_modelinit);
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

    /****************************************************/

#ifdef IWAVE_VERBOSE
    fprintf(stream,"data sampler\n");
    fflush(stream);
#endif

    /* assign default sampler params - still pressure only for this version */
    IASN(tindex,IPNT_0);
    tindex[0]=D_UC;

    if (!((state.model.g.dim > 0) && (state.model.g.dim < RARR_MAX_NDIM+1))) {
      fprintf(stream,"ERROR: main - model dim = %d not in range [1, %d]. ABORT\n",state.model.g.dim,RARR_MAX_NDIM);
      err=E_BADINPUT;
      iwave_destroy(&state);
      abortexit(err,pars,&stream);
    }      
    for (i=0;i<RARR_MAX_NDIM;i++) tmult[i]=REAL_ONE;///((ireal)(state.model.g.dim));
    RASN(scoord,RPNT_0);
  
    /* construct traceterm object */
    strcpy(hdrkey,"hdrfile");
    strcpy(datakey,"datafile");
    err=sampler_construct(&trace,pars,tindex,tmult,scoord,0,hdrkey,datakey,stream);
    if (err) {
      fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
      abortexit(err,pars,&stream);
    }
    /* initial reassignment of source coordinates */
    RASN(scoord,trace.t.tg.src[trace.t.tg.xrec]);

    /* extract trace sampler verbosity flag */
    ps_flint(*pars,"dump_term",&dump_term);

#ifdef IWAVE_VERBOSE
    fprintf(stream,"movie sampler\n");
    fflush(stream);
#endif
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
	fprintf(stream,"ERROR: main from iwave_static_init ",
		"xrec = %d, err = %d. ABORT\n", 
		trace.t.tg.xrec, err);
	abortexit(err,pars,&stream);
      }
      
      /* sampler_init should succeed in record range [first,last) */

      //      fprint_grid(stderr,(state.model).g);
      if (sampler_init(&trace,&(state.model),pars,stream)) {
	fprintf(stream,"ERROR: main from sampler_init. ABORT\n");
	abortexit(err,pars,&stream);
      }

      fprintf(stderr,"IWAVE::asg rkw=%d rk=%d isrc=%d\n",
              retrieveGlobalRank(),retrieveRank(),trace.t.tg.irec);


      if (dump_term) sampler_fprint(&trace, stream);
    
    }

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

  //  quietexit(pars,&stream);

  exit(0);
}

