#include "iwave.h"
/* included to enable access to MAXPATHLEN */
#include <sys/param.h>

/*----------------------------------------------------------------------------*/
/*
  Input file parameter names.
*/
/* static const char PNAMES_OMPNT[] = "omp_nt";       // number of OpenMP threads */
/* static const char PNAMES_NT[] = "nt";              // number of timesteps */
/* static const char PNAMES_MODEL[] = "model";        // model */
/* static const char PNAMES_SRC[] = "srctype";        // source terminator */
/* static const char PNAMES_DUMPPI[] = "dump_pi";     // dump parallel info */
static const char PNAMES_DUMPLDA[] = "dump_lda";   /* dump allocated domain */
static const char PNAMES_DUMPLDC[] = "dump_ldc";   /* dump computational domain */
static const char PNAMES_DUMPLDR[] = "dump_ldr";   /* dump receive domains */
static const char PNAMES_DUMPLDS[] = "dump_lds";   /* dump send domains */
static const char PNAMES_STATS[] = "stats";        /* calculate statistics */
static const char PNAMES_NOPTS[] = "nopts";        /* no post-time-step function */
static const char PNAMES_PRACT[] = "printact";     /* output actions */

/*----------------------------------------------------------------------------*/
/* 
   normal exit - should only be called if everything is OK, so 
   assume stream exists
*/
void quietexit(PARARRAY * pars,FILE ** stream) {
#ifdef IWAVE_VERBOSE
    fprintf(*stream,"Normal exit\n");
#endif
    fflush(*stream);
    fclose(*stream);
    *stream=NULL;
    ps_delete(&pars);
    /*  exit(0); */
}

/*----------------------------------------------------------------------------*/
/* 
   panic - get out quick! clean up environment and exit 
*/
void abortexit(int err,PARARRAY * pars,FILE ** stream) {
    if (*stream) {
	fprintf(*stream,"ABORT on error code %d\n",err);
	fflush(*stream);
	fclose(*stream);
	*stream=NULL;
    }
    else {
	fprintf(stderr,"ABORT on error code %d\n",err);
    }
    ps_delete(&pars);
#ifdef IWAVE_USE_MPI
    MPI_Abort(retrieveGlobalComm(),err);
#else
    exit(1);
#endif
}

/*----------------------------------------------------------------------------*/
/*
  Initializes parallel output stream
*/
int initoutstream(FILE ** stream, int rk, int sz) {

    char *filename; /* workspace for output filename */
    char tmp[20];   /* ditto */
  
    filename = (char *)usermalloc_(MAXPATHLEN * sizeof(char));    
    if ( filename == NULL ) {
	if ( rk == 0 )
	    fprintf(stderr,"IWAVE %dD.\nERROR. Internal: memory allocation for cout filename. ABORT.\n", IWAVE_NDIM);
	return E_BADINPUT;
    }
    if (!getcwd(filename,MAXPATHLEN)) {
      fprintf(stderr,"IWAVE %dD.\nERROR. Internal: getcwd failed for cout filename. ABORT.\n", IWAVE_NDIM);
      userfree_(filename);
      return E_FILE;
    }

    size_t aaa = strlen(filename);
    char bbb = filename[aaa-1];
    if ((aaa > 1) && ( bbb != '/')) 
	strcat(filename,"/");
    if      ( sz < 10 )
	sprintf(tmp,"cout%01d.txt",rk);
    else if ( sz < 100 )
	sprintf(tmp,"cout%02d.txt",rk);
    else if ( sz < 1000 )
	sprintf(tmp,"cout%03d.txt",rk);
    else if ( sz < 10000 )
	sprintf(tmp,"cout%04d.txt",rk);
    else if ( sz < 100000 )
	sprintf(tmp,"cout%05d.txt",rk);
    else
	sprintf(tmp,"cout%d.txt",rk);

    strcat(filename,tmp);

    *stream = fopen(filename, "w");

    if ( *stream == NULL ) {
	if ( rk == 0 )
	    fprintf(stderr,"IWAVE %dD.\nERROR. Internal: cannot open %s. ABORT.\n", IWAVE_NDIM, filename);
	userfree_(filename);
	return E_BADINPUT;
    }

    userfree_(filename);

    /* FROM NOW CAN OUTPUT INTO THIS PROCESSOR'S STREAM */
    storeOutstream(*stream);
    fflush(*stream);
    return 0;
}

/*----------------------------------------------------------------------------*/
/*
  Initialize input parameter structures and populates them.
*/
int readinput(PARARRAY ** pars, FILE * stream, int argc, char **argv) {

    int err=0; /* error code */
    int i;     /* counter */

    /* read parameters */
    if (*pars) {
	fprintf(stderr,"ERROR: iwave::readinput\n");
	fprintf(stderr,"  called on non-null PARARRAY argument\n");
	return E_OTHER;
    }
    *pars = ps_new();
    err=ps_createargs(*pars,argc-1,&(argv[1]));
    if (err) {
	fprintf(stream,
		"ERROR: iwave::readinput from ps_creatargs, err = %d\n",
		err);
	fprintf(stream,"  called with args:\n");
	fprintf(stream,"  argc = %d\n",argc-1);
	for (i=0;i<argc-1;i++) 
	    fprintf(stream,"  argv[%d] = %s\n",i,argv[i+1]);
	return err;
    }
#ifdef IWAVE_VERBOSE
    ps_printall(*(*pars),stderr);
#endif

    return err;
}

/*----------------------------------------------------------------------------*/
/* 
   read miscellaneous flags from param array
*/
void readparpars(PARARRAY * pars,
		 int * stats,
		 int * nopts,
		 int * printact,
		 FILE * out) {

    int err=0;  
    *stats = 1;

    err=ps_flint(*pars, PNAMES_STATS, stats);
    if (err) {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Cannot read %s. Default = %d.\n", PNAMES_STATS, *stats);
#endif
    }
    else {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Found %s = %d.\n", PNAMES_STATS, *stats);
#endif
    }
  
    *nopts = 0;
    if ( ps_flint(*pars, PNAMES_NOPTS, nopts) ) {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Cannot read %s. Default = %d.\n", PNAMES_NOPTS, *nopts);
#endif
    }
    else {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Found %s = %d.\n", PNAMES_NOPTS, *nopts);
#endif
    }
  
    *printact = 0;
    if ( ps_flint(*pars, PNAMES_PRACT, printact) ) {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Cannot read %s. Default = %d.\n", PNAMES_PRACT, *printact);
#endif
    }
    else {
#ifdef IWAVE_VERBOSE
	fprintf(out, "NOTE. Found %s = %d.\n", PNAMES_PRACT, *printact);
#endif
    }
}

/*----------------------------------------------------------------------------*/
/* 
   static storage of MPI comm params 
*/
void storeparallel(PARALLELINFO *pinfo) {
#ifdef IWAVE_USE_MPI
    int ia;
    MPI_Comm_size(pinfo->ccomm, &ia);
    storeSize(ia);
    storeRank(pinfo->lrank);
    storeComm(pinfo->ccomm);
#endif
}

int dump_pi(PARARRAY * par, PARALLELINFO *pinfo, FILE *stream) {
    int flag = 0;
    ps_flint(*par,"dump_pi",&flag);
    if (!flag) return 0;
#ifdef IWAVE_USE_MPI
    int idim, in;/* dimension, neighbor number */
    IPNT offs;   /* neighbor offsets [-1/0/1]  */
    int ndim, nnei;
  
    ndim = pinfo->ndim;
    nnei = pinfo->nnei;
  
    fprintf(stream, "NOTE: Parallel info dump:\n"); 
    fprintf(stream, "NOTE: space dim = %d number of neighbors = %d\n",
	    pinfo->ndim, pinfo->nnei);
    fprintf(stream, "NOTE: Group size = %d, group rank = %d, thread support = %d.\n", 
	    pinfo->wsize, pinfo->wrank, pinfo->threadsupp);
    fprintf(stream, "NOTE: Cartesian grid dims: [%d", pinfo->cdims[0]);
    for ( idim = 1; idim < ndim; ++idim ) fprintf(stream, ",%d", pinfo->cdims[idim]);
    fprintf(stream,"]\n");
    fprintf(stream, "NOTE: this processor Cartesian coordinates: <%d", pinfo->crank[0]);
    for ( idim = 1; idim < ndim; ++idim ) fprintf(stream, ",%d", pinfo->crank[idim]);
    fprintf(stream, ">, rank: %d\nNOTE: send processors [offset]=rk:", pinfo->lrank);
    for ( in = 0; in < nnei; ++in ) {
	if ( gen_i2pnt(ndim, in, offs) ) return E_INTERNAL;
	fprintf(stream, " [%d", -offs[0]); /* -offs, because send offs are subtracted */
	for ( idim = 1; idim < ndim; ++idim ) fprintf(stream, ",%d", -offs[idim]);
	fprintf(stream, "]=%d", pinfo->sranks[in]);
    }
    fprintf(stream, "\nNOTE: recv processors [offset]=rk:");
    for ( in = 0; in < nnei; ++in ) {
	if ( gen_i2pnt(ndim, in, offs) ) return E_INTERNAL;
	fprintf(stream, " [%d", offs[0]);  /* +offs, because recv offs are added */
	for ( idim = 1; idim < ndim; ++idim ) fprintf(stream, ",%d", offs[idim]);
	fprintf(stream, "]=%d", pinfo->rranks[in]);
    }
    fprintf(stream, "\n");
#else
    fprintf(stream, "Parallel info dump: MPI not supported.\n");
#endif
    fflush(stream);
  
    return 0;
}

/*----------------------------------------------------------------------------*/

void dump_ac(PARARRAY * pars, IMODEL * model, FILE * stream) {
    int dump = 0;

    ps_flint(*pars, PNAMES_DUMPLDA, &dump);
    if ( dump ) 
    {
	fprintf(stream, "Allocated ");
	rd_a_dump(&((*model).ld_a), stream);
    }
    dump = 0;
    ps_flint(*pars, PNAMES_DUMPLDC, &dump);
    if ( dump )
    {
	fprintf(stream, "Computational ");
	rd_a_dump(&((*model).ld_c), stream);
    }
    fflush(stream);
}
/*----------------------------------------------------------------------------*/

void dump_rs(PARARRAY * pars, IMODEL * model, FILE * stream, int sends) {
    int dump = 0, nnei, ndim, i, j;
    IPNT ip;
    RDOM *rd;
    //    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    void * fdm = model->specs;
  
    if ( sends ) ps_flint(*pars, PNAMES_DUMPLDS, &dump);
    else         ps_flint(*pars, PNAMES_DUMPLDR, &dump);
  
    if ( dump  && fdm ) {
	if ( sends ) {
	    fprintf(stream, "Send domains dump:\n");
	    rd = (*model).ld_s;
	}
	else {
	    fprintf(stream, "Recv domains dump:\n");
	    rd = (*model).ld_r;
	}
    
	nnei = (*model).nnei;
	ndim = (*model).g.dim;
    
	for ( i = 0; i < nnei; ++i ) {
	    if ( gen_i2pnt(ndim, i, ip) ) continue;
      
	    fprintf(stream, "Neighbor [%d", ip[0]);
	    for ( j = 1; j < ndim; ++j ) fprintf(stream, ",%d", ip[j]);
	    fprintf(stream, "]\n");
      
	    for ( j = 0; j < RDOM_MAX_NARR; ++j ) {
		/*	if (isdyn(fdm,j)) { */
		fprintf(stream, "%d ", j);
		rd_dump(rd + i, j, stream);
		/*	} */
	    }
	    fprintf(stream, "-------------------\n");
	}
    }
    fflush(stream);
}
/*----------------------------------------------------------------------------*/
/*
  Prepare data exchanges - commit MPI data types (MPI_Type_vector)
*/
#ifdef IWAVE_USE_MPI
void prepexch(PARALLELINFO * pinfo, IMODEL * model, IWaveInfo const & ic) {
    /* REMOVED TEMP 13.03.12 WWS
    //,
    //	      int stats, 
    //	      double * wt0) { */
  
    int nnei = (*model).nnei;
    int i, ia;
      
    for ( i = 0; i < nnei; ++i ) {
	for ( ia = 0; ia < RDOM_MAX_NARR; ++ia ) {
	  if (fd_isdyn(ia,ic)) {   
		if ( (pinfo->sranks[i] != MPI_PROC_NULL) && 
		     (pinfo->seinfo[ia][i].buf != NULL) ) 
		    MPI_Type_commit(&(pinfo->seinfo[ia][i].type));
		else pinfo->seinfo[ia][i].type = MPI_DATATYPE_NULL; 
		/* to have valid type */
		if ( (pinfo->rranks[i] != MPI_PROC_NULL) && 
		     (pinfo->reinfo[ia][i].buf != NULL) ) 
		    MPI_Type_commit(&(pinfo->reinfo[ia][i].type));
		else pinfo->reinfo[ia][i].type = MPI_DATATYPE_NULL; 
		/* to have valid type */
	    }
	}
    }
    /*  if ( stats > 0 ) *wt0 = MPI_Wtime(); */
}
#endif

/*-------------------------------------------------------------------------*/
/* REMOVED TEMP - WWS 13.03.12 */
/* prepare timing data at start of time loop */
/*
  #ifdef IWAVE_USE_MPI
  void preptime(double stat_times[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
  int stat_calls[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
  IWAVE * state
  ) {
  int ia, i;
  if ( state->stats > 0 ) state->wt0 = MPI_Wtime();
  
  for ( ia = 0; ia < RDOM_MAX_NARR; ++ia )
  for ( i = 0; i < MODEL_MAX_NACTION; ++i ) {
  stat_times[ia][i][0] = REAL_MAX;
  stat_times[ia][i][1] = -1.0;
  stat_times[ia][i][2] = 0.0;
  stat_calls[ia][i][0] = stat_calls[ia][i][1] = stat_calls[ia][i][2] = 0;        
  }
  }
  #endif
*/
/*-------------------------------------------------------------------------*/
/* REMOVED TEMP - WWS 13.03.12 */
/* output timing data after the end of time loop */
/*
  void printtime(double stat_times[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
  int stat_calls[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
  IWAVE * state, FILE * stream) 
  {

  #ifdef IWAVE_USE_MPI

  double wt1, wtmin, wtmax, actsumtot, actsum;
  int iv, iact, na, ia;
  
  na = state->model.tsinfo.narr;
  
  if ( state->stats > 0 ) fprintf(stream, "NOTE. Statistics.\n");

  if ( state->stats > 0 ) 
  {
  wt1 = MPI_Wtime() - state->wt0;
  MPI_Allreduce(&wt1, &wtmax, 1, MPI_DOUBLE, MPI_MAX, state->pinfo.ccomm);
  MPI_Allreduce(&wt1, &wtmin, 1, MPI_DOUBLE, MPI_MIN, state->pinfo.ccomm);
  fprintf(stream, "    Wall time: %12.4e.\n", wt1);
  fprintf(stream, "    Max  time: %12.4e.\n", wtmax);
  fprintf(stream, "    Min  time: %12.4e.\n", wtmin);
  }
  if ( state->stats > 1 )
  {
  for ( iv = 0; iv < na; ++iv )
  {
  ia = (state->model).tsinfo.arrs[iv];
  fprintf(stream, "    Array: %2d.\n", ia);
            
  for ( iact = 0; iact < MODEL_MAX_NACTION; ++iact )
  {
  if ( stat_calls[ia][iact][2] > 0 )
  {
  fprintf(stream, "      %2d action %s:", ia, MODEL_ACTION_NAMES[iact]);
  fprintf(stream, "  %6d calls,", stat_calls[ia][iact][2]);
  fprintf(stream, "  min = %9.2e (%6d),  max = %9.2e (%6d),  ave = %9.2e\n", 
  (stat_times)[ia][iact][0], stat_calls[ia][iact][0], 
  (stat_times)[ia][iact][1], stat_calls[ia][iact][1], 
  (stat_times)[ia][iact][2] / (double)(stat_calls[ia][iact][2]));
  }
  }
  }

  actsumtot = 0.0;
  for ( iact = 0; iact < MODEL_MAX_NACTION; ++iact )
  {
  actsum = 0.0;
  for ( iv = 0; iv < na; ++iv )
  {
  ia = (state->model).tsinfo.arrs[iv];
  if ( stat_calls[ia][iact][2] > 0 ) actsum += (stat_times)[ia][iact][2] / (double)(stat_calls[ia][iact][2]);
  }
  actsumtot += actsum;
  }
       
  fprintf(stream, "    Action total:         ave = %9.2e, (100%%)\n", actsumtot);
  if ( actsumtot > 0.0 )
  {
  for ( iact = 0; iact < MODEL_MAX_NACTION; ++iact )
  {
  actsum = 0.0;
  for ( iv = 0; iv < na; ++iv )
  {
  ia = (state->model).tsinfo.arrs[iv];
  if ( stat_calls[ia][iact][2] > 0 ) actsum += (stat_times)[ia][iact][2] / (double)(stat_calls[ia][iact][2]);
  }
  fprintf(stream, "      action %s:  ave = %9.2e, (%3d%%)\n", MODEL_ACTION_NAMES[iact], actsum, (int)(actsum / actsumtot * 100.0) );
  }
  }        
  }

  #else
  fprintf(stream,"no timing data in serial mode\n");
  #endif
  }

*/
