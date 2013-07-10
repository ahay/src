#include "sampler.h"

int sampler_construct(SAMPLER * s,
		      PARARRAY * par,
		      const IPNT sindex,
		      const RPNT mindex,
		      const RPNT scoord,
		      int load,
		      const char * hdrkey,
		      const char * datakey,
		      FILE * stream) {
    int err;

    IASN(s->sindex,sindex);
    RASN(s->mindex,mindex);
    RASN(s->scoord,scoord);
    s->load = load;

    /*  fprintf(stderr,"CALL SAMPLER_CONSTRUCT hdrkey=%s datakey=%s\n",hdrkey,datakey);*/
    /* construct traceterm using first comp of sampling and multiplier index
       arrays - these are updated as needed in the run method */
    /* note - sindex, mindex removed from arg list mod of 03.12 */
    err=traceterm_construct(&(s->t),par,
			    /* sindex[0],mindex[0], */
			    load,hdrkey,datakey,stream);
    if (err) {
      fprintf(stream,"Error: sampler_construct from traceterem_construct, err=%d\n",err);
      fprintf(stream,"  hdrkey=%s datakey=%s load=%d\n",hdrkey,datakey,load);
      fprintf(stream,"  PARARRAY:\n");
      ps_printall(*par,stream);
    }
    /*  fprintf(stderr,"EXIT SAMPLER_CONSTRUCT err=%d\n",err);*/
    return err;
}  

int sampler_init(SAMPLER * s,
		 IMODEL * m,
		 PARARRAY * par, 
		 FILE * stream) {
    int err=0;
    err=traceterm_init(&(s->t),m,par,stream);
    if (err) {
	fprintf(stream,"Error: sampler_init from traceterm_init\n");
	fflush(stream);
    }  
    return err;
}

int sampler_destroy(SAMPLER * s) {
    return traceterm_destroy(&(s->t));
}

int sampler_run(SAMPLER * s, IMODEL * m) {
    int i,j;         /* loop counters */
    int err=0;       /* error flag    */
    int rep=1;       /* repeat flag   */

    /* major change 03.12: do same thing on loads and saves */
    /*  if (s->load) { */
    /* load - loop over sindex, also pass corr. multiplier
       both done by modifying trace_term data members - should 
       this be relegated to a trace_term member function?
    */
    for (i=0;i<RARR_MAX_NDIM;i++) {
	if (!err) {
	    /* check that no repeats */
	    rep=1;
	    for (j=0;j<i;j++) 
		if ((s->sindex)[i]==(s->sindex)[j]) rep=0;
	    if (rep) {
		(s->t).index=(s->sindex)[i];
		(s->t).mult =(s->mindex)[i];
		err = traceterm_run(&(s->t),m);
	    }
	}
    }
    /*}
    //else {
    //   save - save only from first field 
    //  (s->t).index=(s->sindex)[0];
    //  err = traceterm_run(&(s->t),m);
    //} */
    
    return err;
}

void sampler_fprint(SAMPLER * s, FILE * fp) {
    int i;
    fprintf(fp,"/*---------------------------------------------------------*/\n");
    fprintf(fp,"IWAVE SAMPLER - uses TRACE HANDLER:\n");
    traceterm_fprint(&(s->t),fp);
    fprintf(fp,"  source, multiplier index arrays:\n");
    for (i=0;i<RARR_MAX_NDIM;i++) 
	fprintf(fp,"  dim=%d sindex=%d mindex=%e\n",i,s->sindex[i],s->mindex[i]);
	  
}

