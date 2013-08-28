/* driver for iwave with elastic staggered grid model
   esg.c
*/
// Modified by Mario Bencomo 2013
// Modifications include:
// 	> new esg time step interfaces
//	> multiple trace output capabilities using MULTI_SAMPLER struct


/* sample of par file for esg */
/*======================================================================

INPUT DATA FOR iwave

------------------------------------------------------------------------
FD:

         order = 2
           cfl = 0.4
          cmin = 1.0
          cmax = 4.5
	  dmax = 5000
	  dmin = 500
         fpeak = 0.010		peak frequency

	 ansol = HI		Analytical solution flag (homogenuous isotropic),
				Requires homogeneous coefficients, point source, 
				and par given refdist.

------------------------------------------------------------------------
Model info:

        lambda = test_lambda.rsf
            mu = test_mu.rsf
       density = test_dn.rsf


------------------------------------------------------------------------
Source info:

      srctype = point/array	
 waveletphase = causal/zerophase
      sampord = 1             	sampling order
       refamp = 1 
      refdist = 1		Required in analytical HI solution case.
        srcin = body		Source injection type. If not included will
				inject in pressure fields.

------------------------------------------------------------------------
Trace info:

      hdrfile = hdr_middle.su
    datafile1 = p0 		z pressure
    datafile2 = p1		x pressure
    datafile3 = v0		z velocity
    datafile4 = s0		zx stress 

    If ansol flag is set, then analytical solution traces are outputted
    for each trace of numerical solution.

------------------------------------------------------------------------
Movie info:

    moviestep = 10.0         time between movie frames (<= 0.0, 
                             or do not include, for no movie)
       movie1 = p0
       movie2 = p1
       movie3 = ansol_p0
       movie4 = ansol_p1

    User must specify the movies outputted for the analytical solutions
    seperate from the numerical solution movies. 

------------------------------------------------------------------------
PML info:

          nl1 = 0.5         z - neg
          nr1 = 0.5         z - pos
          nl2 = 0.5         x - neg
          nr2 = 0.5         x - pos

MPI info:

      mpi_np1 = 1           n_doms along axis 1
      mpi_np2 = 1           n_doms along axis 2
      mpi_np3 = 1           n_doms along axis 3
      partask = 1

Output info:

     printact = 3
      dump_pi = 1           dump parallel/dom. decomp info
     dump_lda = 1           dump grid data for allocated arrays
     dump_ldc = 1           dump grid data for computational arrays
    dump_term = 1           dump terminator data

=========================================================================*/

/* NOTE: There are two versions of esg, Xin Wang's fixed version (with single 
 	 array indexing) and Mario Bencomo's version (width multi index arrays, 
	 and different free surface implementaion ).
	 To use XW's version, uncomment '#define __ESG_STEPS__' and comment out 
	 '#define __ESG_STEPS_MOD__' in esgn.h. For MB's version, do the oppposite. 
	 Another difference is that XW's version should still work when running 
	 in multiple cores, when not including the analytical solution flag 
	 (not sure if it can run in parallel). MB's version still needs to be 
	 checked for parallel case.
*/


/*============================================================================
 *                            BEGIN INCLUDES 
 * ============================================================================*/

#include <iwave.h>
#include <trace_term.h>
#include <pointsrc.h>
#include <sampler.h>
#include <parser.h>
#include <esg_selfdoc.h>
#include <par.h>

#include <esgn.h>
#include <ansol_esgn.h>

#include <esg_movie.h>
#include <ansol_esg_movie.h>

#include <esg_multi_sampler.h> 
#include <ansol_esg_multi_sampler.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NSTR 128
#define IWAVE_VERBOSE
// #define INJ_TIME 250
// #define INJ_TIME_START 260
// #define EARLY_EXIT(x) 2

// #define IWAVE_EXTEND_MODEL

/* uncomment to write to the rk-dep output stream at every time step */
// #define VERBOSE_STEP

/*============================================================================
 *                             END INCLUDES 
 * ============================================================================*/

int main(int argc, char ** argv) {

	int err=0;              /* error flag         */
  	FILE * stream;          /* output stream      */
  	PARARRAY * pars=NULL;   /* parameter array    */
  	IWAVE state;            /* model state        */
	SAMPLER trace;          /* trace sampler      */
	POINTSRC * ptsrc=NULL;  /* acoustic pt src    */
	SAMPLER * arrsrc=NULL;  /* source array       */
	int towed=0;            /* fixed vs. towed    */
	RPNT d;                 /* workspace for i/o  */
	RPNT o;                 /* workspace for i/o  */
	MOVIE mt;               /* movie writer       */
	char * srctype;         /* string length      */
	char * dataval; 	/* string             */
	char hdrkey[NSTR];    	/* key = hdrfile      */
	char datakey[NSTR];    	/* key = datafile     */
	char srckey[NSTR];     	/* key = source       */
	IPNT sindex;           	/* source array index */
	RPNT mult  ;           	/* multiplier array   */
	RPNT scoord;           	/* source cell loc    */
	int dump_term=0;       	/* trace info dump    */
	int istart=0;          	/* start index        */
	int ts;                	/* thread support lvl */
	int rk;                	/* process rank       */

	MULTI_SAMPLER multi_trace;	/* multi-trace sampler    */ //MB
	int multiflag;			/* flag for multi-sampler */ //MB
	
  	PARARRAY      * pars_sol=NULL;
	IWAVE           state_sol;
	SAMPLER         trace_sol;
	MULTI_SAMPLER   multi_trace_sol;
	MOVIE		mt_sol;
	char          * solval;
	char            solkey[NSTR];
	int             solflag;

#ifdef INJ_TIME
	FD_MODEL * fdm;
	FD_MODEL * sol;
	ESGN_TS_PARS * esgnpars;
	ANSOL_ESG_PARS * ansolpars;
#endif

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
#endif
    		fprintf(stderr,"Global MPI_Comm_size = %d\n",retrieveGlobalSize());
#ifdef IWAVE_USE_MPI
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
	//MB	modified
  	err=readinput(&pars,stream,argc,argv);
  	if (err) {
    		fprintf(stderr,"ERROR: main from readinput. ABORT\n");
    		abortexit(err,pars,&stream);
  	}
#ifdef IWAVE_VERBOSE
  	fprintf(stream,"paramtable:\n");
  	fprintf(stream,"--------------------\n");
	ps_printall(*pars,stream);
	fprintf(stream,"--------------------\n");
	fflush(stream);
#endif

#ifdef IWAVE_VERBOSE
  	fprintf(stream,"initparallel_local \n");
  	fflush(stream);
#endif

  	initparallel_local(*pars,stream);

	/******************************************************************
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

	//MB deleted: xargc=argc;
	/* assign key strings */
	strcpy(hdrkey,"hdrfile");
	strcpy(datakey,"datafile");
	strcpy(srckey,"source");
	strcpy(solkey,"ansol");
	
	/* extract trace sampler verbosity flag */
	ps_ffint(*pars,"dump_term",&dump_term);
	
	/* extract source key - default = "point" */
	if (ps_ffcstring(*pars,"srctype",&srctype)) {
		srctype=(char *)usermalloc_(6*sizeof(char)); 
		strcpy(srctype,"point");
	}

	solflag = 0;
	/* extract analytical solution flag */
	if( !ps_ffcstring(*pars,solkey,&solval) ){
		solflag = 1;
		fprintf(stream,"analytical solution flag found: %s = %s\n",solkey,solval);

		//Creating copy of pars for ansol, i.e., pars_sol. Will be needed for movies and traces.
	  	readinput(&pars_sol,stream,argc,argv);
	}



	/************* construct iwave objects *******************************************/
#ifdef IWAVE_VERBOSE
    	fprintf(stream,"iwave_construct\n");
    	fflush(stream);
#endif

    	err=iwave_construct(&state,pars,stream,&esg_modelinit);
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


	/* solflag */
	if (solflag){
#ifdef IWAVE_VERBOSE
    		fprintf(stream,"iwave_construct for ansol\n");
    		fflush(stream);
#endif

		if (!strcmp(solval,"HI")){
			fprintf(stream,"iwave_construct for HI case\n");
    			err=iwave_construct(&state_sol,pars_sol,stream,&ansol_HI_esg_modelinit);
		}
		else { solflag = 0; }
    		if (err) {
      			fprintf(stream,"ERROR: main from iwave_construct for ansol. ABORT\n");
      			iwave_destroy(&state_sol);
      			abortexit(err,pars_sol,&stream);
    		}

#ifdef IWAVE_VERBOSE
    		fprintf(stream,"post-iwave_construct step\n");
    		fflush(stream);
#endif
		err=ansol_link_mdl( stream, &(state_sol.model), &(state.model) );
		if (err) {
      			fprintf(stream,"ERROR: main from ansol_link_mdl. ABORT\n");
      			iwave_destroy(&state_sol);
      			abortexit(err,pars_sol,&stream);
    		}


#ifdef IWAVE_VERBOSE
    		fprintf(stream,"iwave_printf for ansol\n");
    		fflush(stream);
#endif	
    		iwave_printf(&state_sol,pars_sol,stream);

	}/* end solflag */


//printing grid info
#ifdef IWAVE_VERBOSE
	fprintf(stream,"printing grid info:\n");
	fprintf(stream,"for grid g\n");
	fprintf(stream,"     dim = %d\n",state.model.g.dim);
	fprintf(stream,"    gdim = %d\n",state.model.g.gdim);
	int i;
	for(i=0; i< state.model.g.dim; i++){
		fprintf(stream,"axis-%d\n",i);
                fprintf(stream,"         n = %zu\n",state.model.g.axes[i].n);
		fprintf(stream,"         d = %f\n",state.model.g.axes[i].d);
		fprintf(stream,"         o = %f\n",state.model.g.axes[i].o);
                fprintf(stream,"         id = %d\n",state.model.g.axes[i].id);
	}
	fprintf(stream,"for grid gl\n");
	fprintf(stream,"     dim = %d\n",state.model.gl.dim);
	fprintf(stream,"    gdim = %d\n",state.model.gl.gdim);
	for(i=0; i< state.model.g.dim; i++){
		fprintf(stream,"axis-%d\n",i);
                fprintf(stream,"         n = %zu\n",state.model.gl.axes[i].n);
		fprintf(stream,"         d = %f\n",state.model.gl.axes[i].d);
		fprintf(stream,"         o = %f\n",state.model.gl.axes[i].o);
                fprintf(stream,"         id = %d\n",state.model.gl.axes[i].id);
	}
#endif



	/************* construct sampler objects *******************************************/

#ifdef IWAVE_VERBOSE
    	fprintf(stream,"data sampler\n");
    	fflush(stream);
#endif

    	//Looking for datafile
    	multiflag = ps_ffcstring(*pars,"datafile",&dataval);	
	fprintf(stream,"multiflag = %d\n",multiflag);

	/* start multi-trace block */
	if (multiflag){
      	
		//Case where datafile was not found, i.e., multiple output assumed.
#ifdef IWAVE_VERBOSE
    		fprintf(stream,"key \"datafile\" was not found!\n"
                       		"will assume datafile1, datafile2, .... are given\n");
#endif

#ifdef IWAVE_VERBOSE
		fprintf(stream, "pre-construct multi-sampler\n");
#endif
		//esg precons for multi-sampler multi_trace
      		err = esg_multi_sampler_precons( &multi_trace, stream );
    		if (err) {
      			fprintf(stream,"ERROR: main from esg_multi_sampler_precons. ABORT\n");
      			abortexit(err,pars,&stream);
   		}
#ifdef IWAVE_VERBOSE
		fprintf(stream, "construct multi-sampler\n");
#endif
		//constructing multi-sampler multi_trace
		err = multi_sampler_construct( &multi_trace, pars, hdrkey, stream );
		if (err) {
			fprintf(stream,"ERROR: main from multi_sampler_construct. ABORT\n");
			abortexit(err,pars,&stream);
		}
	
		trace = (multi_trace.samplers)[0];

		/* solflag */
		if (solflag){ 
#ifdef IWAVE_VERBOSE
			fprintf(stream, "pre-construct multi-sampler for ansol\n");
#endif
			//ansol_esg precons for multi-sampler multi_trace
      			err = ansol_esg_multi_sampler_precons( &multi_trace_sol, stream );
    			if (err) {
      				fprintf(stream,"ERROR: main from ansol_esg_multi_sampler_precons. ABORT\n");
      				abortexit(err,pars_sol,&stream);
   			}
#ifdef IWAVE_VERBOSE
			fprintf(stream, "construct multi-sampler for ansol\n");
#endif
			//constructing multi-sampler multi_trace
			err = multi_sampler_construct( &multi_trace_sol, pars_sol, hdrkey, stream );
			if (err) {
				fprintf(stream,"ERROR: main from multi_sampler_construct. ABORT\n");
				abortexit(err,pars_sol,&stream);
			}
	
			trace_sol = (multi_trace_sol.samplers)[0];
		}/* end solflag */

	}
	else{ //----- Single trace case, i.e., datafile found ------//

#ifdef IWAVE_VERBOSE
    		fprintf(stream,"key \"datafile\" was found with value: %s\n"
                       		"assigning default sampler params\n",dataval);
#endif

      		/* assign default sampler params */
		IASN(sindex,D_S);
		if (!((state.model.g.dim > 0) && (state.model.g.dim < RARR_MAX_NDIM+1))) {
			fprintf(stream,"ERROR: main - model dim = %d not in range [1, %d]. ABORT\n",state.model.g.dim,RARR_MAX_NDIM);
			err=E_BADINPUT;
			iwave_destroy(&state);
			abortexit(err,pars,&stream);
		}      
		RASN(mult,RPNT_0);
		mult[0]=1;
		RASN(scoord,RPNT_0);

		/* construct traceterm object */
		sindex[0] = D_P0;
		err=sampler_construct(&trace,pars,sindex,mult,scoord,0,hdrkey,datakey,stream);
		if (err) {
			fprintf(stream,"ERROR: main from sampler_construct. ABORT\n");
			abortexit(err,pars,&stream);
		}

		/* solflag */
		if (solflag){ 
			err=sampler_construct(&trace_sol,pars_sol,sindex,mult,scoord,0,hdrkey,datakey,stream);
			if (err) {
				fprintf(stream,"ERROR: main from sampler_construct for ansol. ABORT\n");
				abortexit(err,pars_sol,&stream);
			}
		} /* end solflag */
	}
	/* end multi-trace block */

	/* initial reassignment of source coordinates */
	RASN(scoord,trace.t.tg.src[trace.t.tg.xrec]);



	/************* construct movie objects *******************************************/

#ifdef IWAVE_VERBOSE
    	fprintf(stream,"movie sampler\n");
    	fflush(stream);
#endif
    	err = esg_movie_construct(&mt, stream);
    	if (err) {
      		fprintf(stream,"ERROR: main from movie_construct. ABORT\n");
      		abortexit(err,pars_sol,&stream);
   	}

	/* postconditon: tracegeom trace.tg now equipped with 
	* - number of shot records
	* - shot coordinates for each record
	* - file offset for each record
	*/
	
	/* NOTE: pointsrc object has trivial default construction */

	/* solflag */
	if (solflag){ 
#ifdef IWAVE_VERBOSE
    		fprintf(stream,"movie sampler for ansol\n");
    		fflush(stream);
#endif
    		/* construct movie object */
    		err = ansol_esg_movie_construct(&mt_sol, stream);
    		if (err) {
      			fprintf(stream,"ERROR: main from ansol_esg_movie_construct. ABORT\n");
      			abortexit(err,pars_sol,&stream);
   		}
	} /* end solflag */


	/*****************************************************************************/
	/************* start of trace loop *******************************************/
	/*****************************************************************************/

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

		/************* apply iwave_static_init *****************/
#ifdef IWAVE_VERBOSE
			fprintf(stream,"iwave_static_init\n");
			fflush(stream);
#endif	
		/* NOTE: model coefficients are set here! */
      		/* iwave_static_init should succeed in record range [0,nrec) */
      		err = iwave_static_init( &state,pars,stream,trace.t.tg.xrec,trace.t.tg.first );
      		if(err){  
			fprintf(stream,"ERROR: main from iwave_static_init, xrec = %d, err = %d. ABORT\n", trace.t.tg.xrec, err);
			abortexit(err,pars,&stream);
      		}

		/* solflag */
		if (solflag){ 
#ifdef IWAVE_VERBOSE
			fprintf(stream,"iwave_static_init for ansol\n");
			fflush(stream);
#endif	
      			err = iwave_static_init( &state_sol,pars_sol,stream,trace_sol.t.tg.xrec,trace_sol.t.tg.first );
      			if(err){  
				fprintf(stream,"ERROR: main for ansol from iwave_static_init, xrec = %d, err = %d. ABORT\n", trace_sol.t.tg.xrec, err);
				abortexit(err,pars_sol,&stream);
      			}
#ifdef IWAVE_VERBOSE
			fprintf(stream,"setting medium info for ansol\n");
			fflush(stream);
#endif			

			//setting ansol medium information, after iwave_static_init call!
			err=ansol_HI_esg_medinfo( stream, &(state_sol.model) );
			if (err) {
				fprintf(stream,"ERROR: main from ansol_HI_esg_medinfo. ABORT\n");
				iwave_destroy(&state_sol);
				abortexit(err,pars_sol,&stream);
			}
		} /* end solflag */


		/************* initialize sampler objects *****************/
#ifdef IWAVE_VERBOSE
		fprintf(stream,"initializing sampler objects\n");
		fflush(stream);
#endif	
		/* multi_sampler_init/sampler_init should succeed in record range [first,last) */
		if (multiflag){
      			err = multi_sampler_init( &multi_trace, &(state.model), pars, stream );

		}
		else{
			err = sampler_init( &trace, &(state.model), pars, stream );
		}
		
		if (err){
			fprintf(stream,"ERROR: main from multi_sampler_init/sampler_init. ABORT\n");
			abortexit(err,pars_sol,&stream);
      		}

		if (multiflag) trace = (multi_trace.samplers)[0];
	
		/* solflag */
		if (solflag){ 
			if (multiflag){
				err = multi_sampler_init( &multi_trace_sol, &(state_sol.model), pars_sol, stream );
	
			}
			else{
				err = sampler_init( &trace_sol, &(state_sol.model), pars_sol, stream );
			}
			
			if (err){
				fprintf(stream,"ERROR: main for ansol from multi_sampler_init/sampler_init. ABORT\n");
				abortexit(err,pars_sol,&stream);
			}
	
			if (multiflag){
				trace_sol = (multi_trace_sol.samplers)[0];
			}
		} /* end solflag */ 


      		fprintf(stderr,"IWAVE::esg rkw=%d rk=%d isrc=%d\n",retrieveGlobalRank(),retrieveRank(),trace.t.tg.irec);

      		if (dump_term){
			if (multiflag){
				multi_sampler_fprint( &multi_trace, stream ); 
			}
			else{
 				sampler_fprint(&trace, stream); 
			}
		}


		/************* initialize source *******************************************/
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
			ps_ffint(*pars,"towed",&towed);
			if (!towed) {
				RASN(scoord,RPNT_0);
			}
			/* NOTE THIS PROBABLY ISN'T ANY INTERESTING SOURCE */
			IASN(sindex,D_P);
			err=sampler_construct(arrsrc,pars,sindex,mult,scoord,1,NULL,srckey,stream);
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

		//setting srcinfo for ansol, after samplers, source, and ansol_HI_esg_medinfo are initialized
		/* solflag */
		if (solflag){
#ifdef IWAVE_VERBOSE
			fprintf(stream,"setting srcinfo on ansol\n");
			fflush(stream);
#endif	
			err=ansol_HI_esg_srcinfo( stream, pars_sol, &(state_sol.model), &(trace.t.tg), ptsrc );
			if (err) {
				fprintf(stream,"ERROR: main from ansol_HI_esg_srcinfo. ABORT\n");
				iwave_destroy(&state_sol);
				abortexit(err,pars_sol,&stream);
			}
		} /* end solflag */

		/************* initialize movie *******************************************/
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

		/* solflag */
		if (solflag){ 
#ifdef IWAVE_VERBOSE
      			fprintf(stream,"initialize movie for ansol\n");
      			fflush(stream);
#endif
      			err=movie_init(&mt_sol, &(state_sol.model), pars_sol, &(trace_sol.t.tg), stream);
      			if (err) {
				fprintf(stream,"ERROR: main from movie_init for ansol. ABORT\n");
				abortexit(err,pars_sol,&stream);
      			}
      			if (dump_term) movie_fprint(&mt_sol, stream);
		} /* solflag */



		/************* initialize dynamic fields *******************************************/
#ifdef IWAVE_VERBOSE
      		fprintf(stream,"initialize state\n");
      		fflush(stream);
#endif

		iwave_dynamic_init(&state,istart);
		esg_dynamic_init(&(state.model));

		/* solflag */
		if (solflag){ 
#ifdef IWAVE_VERBOSE
      			fprintf(stream,"initialize state for ansol\n");
      			fflush(stream);
#endif
			/* initialize dynamic fields */
			iwave_dynamic_init(&state_sol,istart);
		} /* end solflag */


#ifdef IWAVE_VERBOSE
		if (solflag) ansol_esg_fprintf( stream, &(state_sol.model) );
#endif



#ifdef INJ_TIME
				fdm = (FD_MODEL *)(state.model.specs);
				sol = (FD_MODEL *)(state_sol.model.specs);

				//Extracting parameters
				esgnpars = (ESGN_TS_PARS *)(fdm->fdpars);
				ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
#endif

		/*************************************************************************/
		/************* start of time loop ****************************************/
		/*************************************************************************/
#ifdef IWAVE_VERBOSE
		fprintf(stream,"time loop\n"); 
		fflush(stream);
#endif

		/* time stepping */


#ifdef EARLY_EXIT
		int mario_count = 0;
		while (state.model.tsind.it <= EARLY_EXIT(trace.t.istop)) {

			fprintf(stderr,"INSIDE: time stepping, it = %d\n",state.model.tsind.it);

// 			if(mario_count>= EARLY_EXIT(trace.t.istop) ){ 
// 				fprintf(stderr,"EARLY EXIT: mario_count = %d\n",mario_count);
// 				exit(1);
// 			}
#endif
#ifndef EARLY_EXIT
		while (state.model.tsind.it <= trace.t.istop) {
#endif

#ifdef INJ_TIME
			if ( state.model.tsind.it<=INJ_TIME ) {
				fdm->tsf = sol->tsf;
				fdm->fdpars = (void *)(ansolpars);
			}
			else if (state.model.tsind.it==INJ_TIME+1) {
				fdm->tsf = esg_step;
				fdm->fdpars = (void *)(esgnpars);
				fprintf(stderr,"after INJ_TIME+1 if\n");
			}
#endif

#ifdef VERBOSE_STEP
			fprintf(stream,"step - update dynamical fields\n");
			fflush(stream);
#endif
			err=iwave_run(&state,stream);
			if (err) {
	  			fprintf(stream,"ERROR: main from iwave_run. ABORT\n");
	  			abortexit(err,pars,&stream);
			}

			/* solflag */
			if (solflag){ 
#ifdef VERBOSE_STEP
				fprintf(stream,"step - update dynamical fields for ansol\n");
				fflush(stream);
#endif
				err=iwave_run(&state_sol,stream);
				if (err) {
	  				fprintf(stream,"ERROR: main from iwave_run for ansol. ABORT\n");
	  				abortexit(err,pars_sol,&stream);
				}
			} /* end solflag */

#ifdef VERBOSE_STEP
			fprintf(stream,"poststep - source\n");
			fflush(stream);
#endif

#ifndef INJ_TIME
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
#endif

#ifdef VERBOSE_STEP
			fprintf(stream,"poststep - traces\n");
			fflush(stream);
#endif
			if (multiflag){
				err=multi_sampler_run( &multi_trace, &(state.model) );
				trace = (multi_trace.samplers)[0];
			}
			else{
				err=sampler_run(&trace,&(state.model));
			}
			if (err) {
				fprintf(stream,"ERROR: main from multi_sampler_run/sampler_run. ABORT\n");
				abortexit(err,pars,&stream);
			}

			/* solflag */
			if (solflag){ 
#ifdef VERBOSE_STEP
				fprintf(stream,"poststep - traces\n");
				fflush(stream);
#endif
				if (multiflag){
					err=multi_sampler_run( &multi_trace_sol, &(state_sol.model) );
					trace_sol = (multi_trace_sol.samplers)[0];
				}
				else{
					err=sampler_run(&trace_sol,&(state_sol.model));
				}
				if (err) {
					fprintf(stream,"ERROR: main from multi_sampler_run/sampler_run for ansol. ABORT\n");
					abortexit(err,pars_sol,&stream);
				}
			} /* end solflag */
		
#ifdef VERBOSE_STEP
			fprintf(stream,"poststep - movie\n");
			fflush(stream);
#endif

			err=movie_run(&mt,&(state.model),stream,0);
			if (err) {
	  			fprintf(stream,"ERROR: main from movie_run. ABORT\n");
	  			abortexit(err,pars,&stream);
			}

			/* solflag */
			if (solflag){ 
#ifdef VERBOSE_STEP
				fprintf(stream,"poststep - movie for ansol\n");
				fflush(stream);
#endif
				err=movie_run(&mt_sol,&(state_sol.model),stream,0);
				if (err) {
	  				fprintf(stream,"ERROR: main from movie_run for ansol. ABORT\n");
	  				abortexit(err,pars_sol,&stream);
				}
			} /* end solflag */

#ifdef VERBOSE_STEP
			fprintf(stream,"update time step info\n");
			fflush(stream);
#endif
			next_step(&(state.model.tsind));
			if (solflag) next_step(&(state_sol.model.tsind));

#ifdef EARLY_EXIT
      			mario_count++;
#endif

      		}
		/*************************************************************************/
		/************* end of time loop ******************************************/
		/*************************************************************************/

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
		if (multiflag){
			err=multi_writetraces( &multi_trace, d, o, stream );
		}
		else{
      			err=writetraces( &(trace.t.tg), d, o, stream );
		}
      		if (err) {
			fprintf(stream,"ERROR: main from tracegeom::writetraces/MULTI_SAMPLER::multi_writetraces. ABORT\n");
			abortexit(err,pars,&stream);
		}

		/* solflag */
		if (solflag){ 
#ifdef IWAVE_VERBOSE
      			fprintf(stream,"writetraces for ansol\n");
      			fflush(stream);
#endif
			if (multiflag){
				err=multi_writetraces( &multi_trace_sol, d, o, stream );
			}
			else{
      				err=writetraces( &(trace_sol.t.tg), d, o, stream );
			}
      			if (err) {
				fprintf(stream,"ERROR: main from tracegeom::writetraces/MULTI_SAMPLER::multi_writetraces for ansol. ABORT\n");
				abortexit(err,pars_sol,&stream);
			}
		} /* end solflag */

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
	/*****************************************************************************/
	/************* end of trace loop *********************************************/
	/*****************************************************************************/
		
/*	if (srctype) userfree_(srctype);*/
	
#ifdef IWAVE_VERBOSE
	//printing pml RDOM
	esg_printpml( pars, stream, &(state.model) );
#endif
	
	/************* destroy static objects and exit *****************************************/
#ifdef IWAVE_VERBOSE
    	fprintf(stream,"movie_destroy\n");
    	fflush(stream);
#endif

    	movie_destroy(&mt);
	if (solflag) movie_destroy(&mt_sol);

#ifdef IWAVE_VERBOSE
    	fprintf(stream,"multi_sampler_destroy/sampler_destroy\n");
    	fflush(stream);
#endif
	if (multiflag){
		multi_sampler_destroy( &multi_trace );
		if (solflag) multi_sampler_destroy( &multi_trace_sol );
	}
	else
	{
    		sampler_destroy(&trace);
		if (solflag) sampler_destroy(&trace_sol);

	}
#ifdef IWAVE_VERBOSE
    	fprintf(stream,"iwave_destroy\n");
    	fflush(stream);
#endif

    	iwave_destroy(&state);
	if (solflag) iwave_destroy(&state_sol);
	//MB added:
   	iwave_fdestroy();

#ifdef IWAVE_USE_MPI
	} /* end nontriv comm branch */

#ifdef IWAVE_VERBOSE
  	fprintf(stream,"MPI_Finalize\n");
  	fflush(stream);
#endif /* IWAVE_VERBOSE */

  	MPI_Barrier(MPI_COMM_WORLD);
  	MPI_Finalize();
#endif /* IWAVE_USE_MPI */

#ifdef IWAVE_VERBOSE
  	fprintf(stream,"quietexit\n");
  	fflush(stream);
#endif

  	quietexit(pars,&stream);
  	if (solflag) quietexit(pars_sol,&stream);

  	exit(0);
}

