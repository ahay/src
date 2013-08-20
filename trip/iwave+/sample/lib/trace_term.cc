#include "trace_term.h"

int traceterm_construct(TRACE_TERM * tr, 
			PARARRAY * par,
			int load,
			const char * hdrkey,
			const char * datakey,
			FILE * stream) {
  
  /* workspace for initialization of trace geometry */
  int err=0;
  char * hfile;          /* input data file name - extracted from input */
  char * dfile;          /* input data file name - extracted from input */
  
  /*  fprintf(stderr,"CALL TRACETERM_CONSTRUCT hdrkey=%s datakey=%s\n",hdrkey,datakey);*/
  /* read data file - if no output, read null */
  hfile=NULL;
  dfile=NULL;
  /* 09.12.09: added keys as parameters; data key must always be non-void
     string, but hdr key may not be */
  /* 13.04.13: force hdrkey to be present. This is a big change: no longer
     possible to specify simulation with just one input/output file. Must 
     supply separate header file.
  */
  if (!datakey || !hdrkey) {
    fprintf(stream,"ERROR: traceterm_construct\n");
    fprintf(stream,"  either key hdrfile (source of trace headers)\n");
    fprintf(stream,"  or key datafile (data trace input/output) not provided\n");
    return E_FILE;
  }
  ps_flcstring(*par,hdrkey,&hfile);
  ps_flcstring(*par,datakey,&dfile);

  /* major change for this edition: if just a datafile name is given, then
     the header and data files are presumed to be the same, and the data file
     is used to supply headers. */
  /* major change 13.04.13: BOTH header and data filenames must be present in parfile.*/
  if (!hfile) {
    fprintf(stream,"ERROR: traceterm_construct\n");
    fprintf(stream,"failed to extract header (prototype) filename from partable");
    return E_FILE;
  }
  if (!dfile) {
    fprintf(stream,"ERROR: traceterm_construct\n");
    fprintf(stream,"failed to extract data filename from partable");
    return E_FILE;
  }

  /* construct trace geometry object */
  err=construct_tracegeom(&(tr->tg),hfile,dfile,SRC_TOL,stream);
  
  /*  fprintf(stderr,"return from construct_tracegeom\n");*/
  
  if (hfile) userfree_(hfile);
  if (dfile) userfree_(dfile);
  if (err) {
    fprintf(stream,"ERROR: traceterm_construct\n");
    fprintf(stream,"return from construct_tracegeom with err=%d\n",err);
    return err;
  }
    
  /* read sampling order */
  tr->order=0;
  ps_flint(*par,"sampord",&(tr->order));

  /* record load flag */
  tr->load=load;

  /* record sampling and multipler indices - must be sanity-checked
     in run method, since they may be overwritten by the sampler.
  */
  /* modification of 03.12 - these never made sense here, since they
     are set in the run method and controlled one level up, by
     sampler. So set to default values. */
  /*  tr->index=sindex; */
  /*  tr->mult=mindex; */
  tr->index=0;
  tr->mult=REAL_ONE;

  /* default initialization of istart, istop */
  tr->istart=0;
  tr->istop=0;

  /*  fprintf(stderr,"EXIT TRACETERM_CONSTRUCT err=%d\n",err);  */
  return err;
}

int traceterm_init(TRACE_TERM * tr, 
		   IMODEL * m, 
		   PARARRAY * par, 
		   FILE * stream) {

  int err=0;                    /* error flag */
  int usernt=0;                 /* user nt override */
  float usert0=0.0;             /* user t0 override */
  IPNT n;                       /* axis lengths, local grid */
  RPNT o;                       /* axis origins, local grid */
  RPNT og;                      /* axis origins, global grid */
  RPNT d;                       /* axis steps, local grid */
  IPNT axord;                   /* axis order array */

  /* extract grid params - space and time */
  get_n(n,m->gl);
  get_o(o,m->gl);
  get_d(d,m->gl);
  get_o(og,m->g);
  get_ord(axord,m->g);

  /* user overrides */
  ps_flint(*par,"nt",&usernt);
  ps_flreal(*par,"t0",&usert0);
  
  /* initialize tracegeom */
  err=init_tracegeom(&(tr->tg),
		     og,n,d,o,axord,
		     tr->order,
		     m->tsind.dt,(m->g).dim,
		     usernt,usert0,tr->load,
		     stream);

#ifdef IWAVE_VERBOSE
  fprintf(stream,"traceterm_init: (tr->tg).t0 = %e, (m->tsind).dt = %e\n",(tr->tg).t0,(m->tsind).dt);
  fflush(stream);
#endif
  /* start and stop from delrt, ns - OK even if (err)*/
  tr->istart=(int)((((tr->tg).t0)/((m->tsind).dt))+0.1);
  tr->istop=tr->istart + (tr->tg).nt - 1;

  if (err)
    fprintf(stream,"Error from traceterm_init: err=%d\n",err);
  
  return err;
}

int traceterm_destroy(TRACE_TERM * tr) {
  destroy_tracegeom(&(tr->tg));
  return 0;
}

int traceterm_run(TRACE_TERM * tr, IMODEL * m) {
  /* workspace for array geometry */
  IPNT n, gs, ge, n0, gs0, ge0;

  /* sanity-check sample array index */
  if (tr->index < 0 || tr->index > (m->ld_a).narr-1) {
    return E_BADINPUT;
  }

  /* NO-OP only on step boundary */
  if ((m->tsind).iv != 0) return 0;

  /* read LOCAL grid params from grid to be sampled */
  ra_gse(&((m->ld_a)._s[tr->index]),gs0,ge0);
  ra_size(&((m->ld_a)._s[tr->index]), n0);
  ra_gse(&((m->ld_c)._s[tr->index]),gs,ge);
  ra_size(&((m->ld_c)._s[tr->index]), n);

  /* if start is reached, start recording, incrementing counter */
  if ((m->tsind).it >  tr->istart-1 &&
      (m->tsind).it <=  tr->istop      )

    sampletraces(&(tr->tg),
		 tr->order,
		 tr->load,
		 (m->tsind).it - tr->istart,
		 n0, gs0,
		 n, gs,
		 (m->ld_a)._s[tr->index]._s,
		 tr->mult);

  return 0;
}

void traceterm_fprint(TRACE_TERM const * tr, FILE * fp) {
  fprintf(fp,"/*---------------------------------------------------------*/\n");
  fprintf(fp,"TRACE SAMPLER\n");
  fprintf(fp,"istart   = %d\n",tr->istart);
  fprintf(fp,"istop    = %d\n",tr->istop);
  fprintf(fp,"order    = %d\n",tr->order);
  fprintf(fp,"load     = %d\n",tr->load);
  fprintf(fp,"index    = %d\n",tr->index);
  fprintf(fp,"mult     = %e\n",tr->mult);
  fprint_tracegeom(&(tr->tg),fp);
}
