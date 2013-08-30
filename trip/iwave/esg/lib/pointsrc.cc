#include "pointsrc.h"

#define DT_DBG

// #define IWAVE_VERBOSE

/* computing approximation of 1D delta function at given point x, \
   where delta function is centered at 0 with supp = [-1,1], and total energy = 1 */
/*----------------------------------------------------------------------------*/
ireal approx_delta( int order, ireal x ){
/*----------------------------------------------------------------------------*/
// #ifdef IWAVE_VERBOSE
// 	fprintf(stderr,">>>> inside approx_delta\n");
// 	fprintf(stderr," order = %d, x = %g\n",order,x);
// #endif
	float val = -1;
	if ( order==0 ){
		if (fabs(x)<=1) val = 0.5;
		else val = 0;
	}
	else if ( order==1 ){
		float a0 = 1.0;
		float a1 = -a0;
		if (fabs(x)< 1 ){
			if (x>=0) val = a1*x + a0;
			else val = -a1*x + a0;
		}
		else val = 0;
	}
	else if ( order==2 ){
		float a0 = 2.4;
		float a1 = -a0*1.5;
		float a2 = a0*0.5;
		if (fabs(x)< 1){
			if (x>=0) val = a2*x*x + a1*x + a0;
			else val = a2*x*x - a1*x + a0;
		}
		else val = 0;
	}
// #ifdef IWAVE_VERBOSE
// 	fprintf(stderr,"val=%g\n",val);
// #endif
	return val;
}

/* point sampling */
/*----------------------------------------------------------------------------*/
int pointsource(IPNT is,
		RPNT rs,
		int order,
		ireal scramp,
		ireal src,
		RARR arr,
		int arr_index ) {
/*----------------------------------------------------------------------------*/
// #ifdef IWAVE_VERBOSE
// 	fprintf(stderr,">>> Inside pointsource\n");
// 	fprintf(stderr,"          src = %e\n",src);
// 	fprintf(stderr," src*scramp = %e\n",src*scramp);
// #endif
  	IPNT gs, ge;
  	ireal fac;
  	int ndim;
  	ireal p;

  	FILE *stream;
  	stream = retrieveOutstream();

  	ra_ndim(&arr, &ndim);
  	ra_gse(&arr, gs, ge);

	IPNT loc_is;
	RPNT loc_rs;
	IASN(loc_is,is);
	RASN(loc_rs,rs);

	/* accounting for duality of velocity grids */
	if (arr_index==D_V0){
		if ( fabs(loc_rs[0]-0.5)<1e-9 ) {
	 		loc_rs[0] = 0;
		}
		if ( loc_rs[0]>0.5 ) {
			loc_rs[0] -= 0.5;
		}
		if ( loc_rs[0]<0.5 ){
			loc_rs[0] += 0.5;
			loc_is[0]--;	
		}
	}
	if (arr_index==D_V1){
		if ( fabs(loc_rs[1]-0.5)<1e-9 ) {
	 		loc_rs[1] = 0;
		}
		if ( loc_rs[1]>0.5 ) {
			loc_rs[1] -= 0.5;
		}
		if ( loc_rs[1]<0.5 ){
			loc_rs[1] += 0.5;
			loc_is[1]--;	
		}
	}
	if (arr_index==D_V2){
		if ( fabs(loc_rs[2]-0.5)<1e-9 ) {
	 		loc_rs[2] = 0;
		}
		if ( loc_rs[2]>0.5 ) {
			loc_rs[2] -= 0.5;
		}
		if ( loc_rs[2]<0.5 ){
			loc_rs[2] += 0.5;
			loc_is[2]--;	
		}
	}

// #ifdef IWAVE_VERBOSE
// 	int idim;
// 	fprintf(stderr,"for arr_index = %d\n",arr_index);
// 	for (idim=0; idim<ndim; idim++)
// 		fprintf(stderr,"loc_is[%d] = %d, loc_rs[%d] = %g\n",idim,loc_is[idim],idim,loc_rs[idim]);
// #endif

	if (ndim==1){
		/* [0] */
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 ){
			fac = approx_delta( order, -loc_rs[0] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		/* [1] */
		loc_is[0]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 ){
			fac = approx_delta( order, 1.0-loc_rs[0] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[0]--;
	}
	if (ndim==2) {
		/* [0,0] */
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
/*			fprintf(stderr,"for 0,0: p + src*scramp*fac = %g\n",p + src*scramp*fac);*/
		}

		/* [1,0] */
		loc_is[0]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
// 			fprintf(stderr,"for 1,0: p + src*scramp*fac = %g\n",p + src*scramp*fac);

		}
		loc_is[0]--;

		/* [0,1] */
		loc_is[1]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
/*			fprintf(stderr,"for 0,1: p + src*scramp*fac = %g\n",p + src*scramp*fac);*/
		}
		loc_is[1]--;

		/* [1,1] */
		loc_is[0]++; loc_is[1]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
/*			fprintf(stderr,"for 1,1: p + src*scramp*fac = %g\n",p + src*scramp*fac);*/
		}
		loc_is[0]--; loc_is[1]--;
	}
	if (ndim==3){
		/* [0,0,0] */
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			fac *= approx_delta( order, -loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}

		/* [1,0,0] */
		loc_is[0]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			fac *= approx_delta( order, -loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[0]--;

		/* [0,1,0] */
		loc_is[1]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			fac *= approx_delta( order, -loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[1]--;

		/* [0,0,1] */
		loc_is[2]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			fac *= approx_delta( order, 1.0-loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[2]--;

		/* [1,1,0] */
		loc_is[0]++; loc_is[1]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			fac *= approx_delta( order, -loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[0]--; loc_is[1]--;

		/* [0,1,1] */
		loc_is[1]++; loc_is[2]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, -loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			fac *= approx_delta( order, 1.0-loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[1]--; loc_is[2]--;

		/* [1,0,1] */
		loc_is[0]++; loc_is[2]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, -loc_rs[1] );
			fac *= approx_delta( order, 1.0-loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[0]--; loc_is[2]--;

		/* [1,1,1] */
		loc_is[0]++; loc_is[1]++; loc_is[2]++;
		if ( loc_is[0]>gs[0]-1 && loc_is[0]<ge[0]+1 &&
		     loc_is[1]>gs[1]-1 && loc_is[1]<ge[1]+1 &&
		     loc_is[2]>gs[2]-1 && loc_is[2]<ge[2]+1 ){
			fac  = approx_delta( order, 1.0-loc_rs[0] );
			fac *= approx_delta( order, 1.0-loc_rs[1] );
			fac *= approx_delta( order, 1.0-loc_rs[2] );
			p = ra_gget( &arr, loc_is );
			ra_gset( &arr, loc_is, p + src*scramp*fac );
		}
		loc_is[0]--; loc_is[1]--; loc_is[2]--;
	}
 
  	return 0;
}

/*----------------------------------------------------------------------------*/
int pointsrc_init(POINTSRC * tr, IMODEL * m, PARARRAY * par, tracegeom *tg, FILE * stream) {
/*----------------------------------------------------------------------------*/
	int     err = 0;  	/* error flag */
	char  * wp;          	/* wavelet phase */
	int     iw;           	/* half-width */
	int     i;            	/* counter */
	int     ndim;        	/* dimension of grids */
	RPNT    d;             	/* steps */
	ireal   prod_d;      	/* cell volume scale factor for delta */
	char  * srcfile;    	/* workspace for file name */
	segy    trsrc;    	/* segy workspace for reading file */
	Value   val;      	/* workspace for reading segy headers */
	int     tmpnt;        	/* length of time series read from file */
	int     lnt;         	/* length of extended time series, for integration */
	ireal   tmpdt;    	/* time step for time series read from file */
	ireal   tmpt0;    	/* time origin for time series read from file */
	ireal   tmax;       	/* max time for either trace or source time series */
	int     wl;         	/* length of workspace for cubic spline interp */
	ireal * wk;     	/* workspace for cubic spline interp */
	int     iend = 1;   	/* end condition for cubic spline interp */
	ireal   t0;       	/* time origin on model dt grid */
	ireal   tmp0, tmp1, q; 	/* workspace for in-place trapezoidal rule */
	IPNT    tis;        	/* buffer for coefficient sampling near source */
	int     iflag;        	/* test flag */
	ireal   tdt;      	/* ireal buffer for dt */
	
	ireal refvel;      	/* reference velocity for target prop wavelet */
	ireal refbou;		/* reference bouyancy for target prop wavelet */
	ireal refkappa;        	/* near-source bulk modulus */
	ireal refdist;      	/* reference distance for target prop wavelet */
	ireal refamp;         	/* reference amplitude for target prop wavelet */
	ireal fpeak;         	/* peak frequency for Ricker computation */
	
	segy    trdbg;       	/* workspace for building output segy */
	IPNT    gs, ge;       	/* workspace for bouyancy exchange comps */
	ireal * resc;

// 	int t_grid_off = 1;

	/* extracting srcin flag */
	char *srcin_val;
	tr->srcin_flag = 0;
	if ( !ps_ffcstring(*par,"srcin", &srcin_val) && !strcmp("body",srcin_val) ) 
		tr->srcin_flag = 1;

  	/* MPI workspace */
#ifdef IWAVE_USE_MPI
  	int rk, sz;
  	MPI_Comm cm;
  	ireal *procbuf;
#endif
  
  	/* end declarations */

  	stream = retrieveOutstream();

	/* assign default reference values */
	refvel   = CREF_DEF;
	refbou   = REAL_ZERO;
	refkappa = REAL_ZERO;
	refdist  = RREF_DEF;
	refamp   = REAL_ONE;
	fpeak    = FPEAK_DEF;        
	
	/* get array of grid steps */
	get_d(d, m->gl);
	
	/* extract dimension */
	rd_ndim(&m->ld_a, D_MP00, &ndim);
	
	/* TV */
	IASN(tr->is,IPNT_0);
	RASN(tr->rs,RPNT_0);
	tr->is[0]=tg->is[0]; tr->rs[0]=tg->rs[0]; tis[0]=tr->is[0];
	if (ndim > 1) { tr->is[1]=tg->is[1]; tr->rs[1]=tg->rs[1]; tis[1]=tr->is[1]; }
	if (ndim > 2) { tr->is[2]=tg->is[2]; tr->rs[2]=tg->rs[2]; tis[2]=tr->is[2]; } 
	
	/* read sampling order */
	tr->order = 0;
	ps_ffint(*par, "sampord", &(tr->order));
	
	/* either read from file, or it's a gaussian */
	tr->fpsrc = NULL;
	tr->fpdbg = NULL;
	
	/* extract near-source bulk modulus from grid - necessary for
	   several (but not all) cases, and not enough of an expense to 
	   avoid.
	
	   WWS, 29.11.08: to avoid failure if integer part of source point is in 
	   physical grid (including physical boundary points) but not in comp.
	   grid, check neighboring grid points - in that case, at least one
	   of these should be in part of comp grid, where kappa is initialized.
	*/
	
	rd_gse(&(m->ld_a), D_MP00, gs, ge);
	for ( i = 0; i < ndim; ++i ) {
		if ( (tr->is[i] < gs[i]) && (tr->is[i]+1 > gs[i]-1) ) tis[i]++;
		if ( (tr->is[i] > ge[i]) && (tr->is[i]-1 < ge[i]+1) ) tis[i]--;
	} 
	iflag=1;
	for (i = 0; i < ndim; ++i) 
		if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;
	
	if (iflag) refkappa = rd_gget(&(m->ld_a), D_MP00, tis);
	
	fprintf(stream,"NOTE: proc=%d sample kappa at ",retrieveRank());
	for (i=0;i<ndim;++i) 
		fprintf(stream,"index[%d]=%d ",i,tis[i]);
	fprintf(stream,"sample flag=%d\n",iflag);
	
	/* extract near-source bouyancy from grid - necessary for
	   several (but not all) cases, and not enough of an expense to 
	   avoid.
	
	   WWS, 04.03.09: need to do this here too - to avoid failure if
	   integer part of source point is in physical grid (including
	   physical boundary points) but not in comp.  grid, check
	   neighboring grid points - in that case, at least one of these
	   should be in part of comp grid, where kappa is initialized.
	*/
	
	rd_gse(&(m->ld_a), D_MV0, gs, ge);
	for ( i = 0; i < ndim; ++i ) {
		tis[i] = tr->is[i];
		if ( (tr->is[i] < gs[i]) && (tr->is[i]+1 > gs[i]-1) ) tis[i]++;
		if ( (tr->is[i] > ge[i]) && (tr->is[i]-1 < ge[i]+1) ) tis[i]--;
	} 
	iflag=1;
	for (i = 0; i < ndim; ++i) 
		if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;
	
	if (iflag) refbou = rd_gget(&(m->ld_a), D_MV0, tis);
		fprintf(stream,"NOTE: proc=%d sample bouyancy at ",retrieveRank());
	for (i=0;i<ndim;++i) 
		fprintf(stream,"index[%d]=%d ",i,tis[i]);
	fprintf(stream,"sample flag=%d\n",iflag);
	
#ifdef IWAVE_USE_MPI
	rk = retrieveRank();
	sz = retrieveSize();
	cm = retrieveComm();
	if ( rk == 0 ) {
		procbuf = (ireal*)usermalloc_(sz * sizeof(ireal));
		if ( procbuf == NULL ) return E_ALLOC;
	}
	MPI_Gather(&refkappa, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
	if(rk ==0 ) {
		refkappa = 0.0;
		for ( i = 0; i < sz; ++i ) refkappa=iwave_max(refkappa,procbuf[i]);
	}
	MPI_Bcast(&refkappa, 1, IWAVE_MPI_REAL, 0, cm);
	MPI_Gather(&refbou, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
	if(rk ==0 ) {
		refbou = 0.0;
		for ( i = 0; i < sz; ++i ) refbou=iwave_max(refbou,procbuf[i]);
		userfree_(procbuf);//RN
	}
	MPI_Bcast(&refbou, 1, IWAVE_MPI_REAL, 0, cm);
#endif	
	
	if (refbou > REAL_ZERO) {
		fprintf(stream,"NOTE: in pointsrc, using  bouyancy at source location = %e\n", refbou);
	}
	else {
		fprintf(stream,"ERROR: in pointsrc, ref bouyancy nonpositive, = %e\n",refbou);
		return E_OUTOFBOUNDS;
	}
	
	if (refkappa > REAL_ZERO) {
		fprintf(stream,"NOTE: in pointsrc, using  bulk mod at source location = %e\n", refkappa);
	}
	else {
		fprintf(stream,"ERROR: in pointsrc, ref bulk mod nonpositive, = %e\n",refkappa);
		return E_OUTOFBOUNDS;
	}
	
	/* get reference velocity, either from parameters or from bouyancy
	   and bulk modulus near source location. Since velocity is not
	   stored in RDOM, must be computed from bouyancy and bulk mod at
	   slightly different points - this of course does not matter if
	   source is located in homogeneous region.
	*/
	
	if (ps_ffreal(*par,"refvel",&refvel)) {  
	
		/* compute velocity from bouyancy and bulk modulus */
	
		refvel = sqrt(refkappa * refbou);
	
		fprintf(stream,"NOTE: in pointsrc, using velocity computed from \n");
		fprintf(stream,"      bulk mod and bouyancy near source location; \n");
		fprintf(stream,"      computed value = %e\n", refvel);
	
	}
	else {
		fprintf(stream,"NOTE: in pointsrc, using velocity from param table = %e\n", refvel);
	}
	
	/* Either read reference distance from parameters, or use default. */
	ps_ffreal(*par,"refdist",&refdist);
	if (refdist>0) {
		fprintf(stream,"NOTE: in pointsrc, using reference distance = %e\n", refdist);
	}
	else {
		fprintf(stream,"NOTE: in pointsrc, read nonpos. reference distance = %e\n.", refdist);
// 		fprintf(stream,"      this implies that wavelet will be read from file and used \n");
// 		fprintf(stream,"      directly on RHS as multiplier of spatial delta, rather than\n");
// 		fprintf(stream,"      to produce target propagating pulse.\n");
	}
	
	/* Either read reference amplitude from parameters, or use default. */
	
	ps_ffreal(*par,"refamp",&refamp);
	fprintf(stream,"NOTE: in pointsrc, using reference amplitude = %e\n", refamp);
	
	/* read peak frequency from parameters, or use default (only used in option II) */
	if (ps_ffreal(*par,"fpeak", &fpeak)) {
		fprintf(stream,"NOTE: pointsrc_init - using default ");
		fprintf(stream,"peak frequency (fpeak) = %e\n",fpeak);
	}
	fprintf(stream,"NOTE: pointsrc_init, using peak frequency (fpeak) = %e\n",fpeak);

	/* setting src_d */
	RASN( tr->src_d, RPNT_0);
	if (tr->srcin_flag) {
		//default value
		tr->src_d[0] = REAL_ONE;
	}


	/* Option I: read source from file */
	if (!ps_ffcstring(*par,"source", &srcfile))  {
		if (!(tr->fpsrc = iwave_const_fopen(srcfile, "r",NULL,stream))) {
			fprintf(stream, "Error: pointsrc_init - failed to open source file\n");
			return E_FILE;
		}
		if (fseek(tr->fpsrc,0L,SEEK_SET)) {
			fprintf(stream,"Error: pointsrc_init - failed to seek to start of file\n");
			return E_FILE;
		}
		if (!fgettr(tr->fpsrc, &trsrc)) {
			fprintf(stream,"Error: pointsrc_init - failed to read source file\n");
			return E_FILE;
		}
		iwave_fclose(tr->fpsrc);

		/* at this point, the data member of trsrc contains the source wavelet
		   at an external sample rate - read headers relevant for single dilat
		   point source.
		*/
		gethdval(&trsrc, "ns", &val);
		tmpnt = vtoi(hdtype("ns"), val);
		gethdval(&trsrc, "dt", &val);
		tmpdt = 0.001 * vtof(hdtype("dt"), val);
		gethdval(&trsrc, "delrt", &val);
		tmpt0 = vtof(hdtype("delrt"), val);		
		
		/* calculate istart, length of resamp wavelet,
		   allocate buffer to hold it.
		*/
		tr->istart = (int)(tmpt0/((m->tsind).dt));
		t0 = (m->tsind).dt * tr->istart;
		/*    tr->n = (int)(tmpnt*tmpdt/((m->tsind).dt))+1;*/

		/* final time for ext source set to max (final time for
		   trace, final time for input source) */
		tmax = iwave_max(tmpt0 + tmpnt * tmpdt, tg->t0 + tg->nt * ((m->tsind).dt));
		/* ext src array length set to number of samples from istart
		   to max time */
		tr->n = (int)(tmax/((m->tsind).dt)) + 1;
		tr->w = (ireal *)usermalloc_(sizeof(ireal)*(tr->n));
		
		/* allocate buffer for integrated wavelet 
		   at input sample rate */
		lnt  = (int)(tr->n * ((m->tsind).dt) / tmpdt) + 1;
		resc = (ireal *)usermalloc_(sizeof(ireal) * lnt);
		for (i = 0; i < tmpnt; i++) resc[i] = trsrc.data[i];
		for (i = tmpnt; i < lnt; i++) resc[i] = trsrc.data[tmpnt-1];
		
		/* interpolation workspace */
		wl = cubic_getworksize(lnt);
		wk = (ireal *)usermalloc_(sizeof(ireal) * wl);
		if (!wk) return E_ALLOC;
		
		fprintf(stream,"NOTE: source wavelet is read from file, i.e., option I.\n");
		/* Option Ia: normalize to produce propagating wavelet at distance
		   refdist, speed refvel, if refdist > 0. */
		if (refdist > 0) {
		
			fprintf(stream,"NOTE: in pointsrc_init, compute wavelet for target pulse (Option Ia)\n");
			
			/* integrate once - trapezoidal rule in-place 
				w[j] <- sum_{i=1}^{i=j} 0.5*dt*(w[i-1]+w[i])
			   leave multiplication by dt for scaling step
			*/
			
			tmp0 = resc[0];
			q = 0.0;
			resc[0] = 0.0;
			for (i = 1; i < lnt; i++) {
				tmp1 = resc[i];
				q += 0.5 * tmpdt * (tmp0 + tmp1);
				resc[i] = q;
				tmp0 = tmp1;
			}
			
			tdt = (ireal)((m->tsind).dt);
			/* interpolate */
			if ((err=cubic_(&tmpt0, &tmpdt, resc, &lnt, &t0, &(tdt), tr->w, &(tr->n), &iend,wk,&wl))) {
				fprintf(stream,"Error: pointsrc_init - from cubic\n");
				return err;
			}
		}
		
		/* Option Ib: use wavelet as is, if refdist <= 0 */
		else {
		
			fprintf(stream,"NOTE: in pointsrc_init, using file wavelet directly on RHS (Option Ib)\n");
			tdt = (ireal)((m->tsind).dt);
		
			/* interpolate */
			if ((err=cubic_(&tmpt0, &tmpdt, resc, &lnt, &t0, &(tdt), tr->w,&(tr->n),&iend,wk,&wl))) {
				fprintf(stream,"Error: pointsrc_init - from cubic\n");
				return err;
			}
		}
		
		/* clean up */
		userfree_(wk);//RN
		userfree_(srcfile);//RN
		userfree_(resc);//RN
	}
	
	/* Option II, generate source */
	else {
		if(ps_ffreal(*par,"refdist",&refdist)){
			/* check that reference distance is positive - only legal option */
			if (!(refdist>0)) {
				fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
				fprintf(stream,"negative reference distance = %e no legal in this case\n",refdist);
				return E_OTHER;
			}
		}
		else{
			/* if refdist not found, set to default */
			refdist = RREF_DEF;
			fprintf(stream,"Error: refdist not found, using default value!\n");
		} 
		fprintf(stream,"NOTE: in pointsrc_init, computing Gaussian derivative RHS for target Ricker pulse at f=%e, r=%e\n",fpeak,refdist);
	
		/* Option II(a): create Ricker wavelet for body force injection of source
		with specified peak frequency, amplitude at spec'd distance. */
		if (tr->srcin_flag) {
			fprintf(stream,"NOTE: injecting source as a body force, i.e., source present in vel. updates.\n");
			/* RHS wavelet is Ricker, i.e., second derivative of gaussian */
			tr->w = getrick( &iw, (m->tsind).dt, fpeak );
			tr->n = 2 * iw+1;
#ifdef IWAVE_VERBOSE
fprintf(stderr,">>>> Inside pointsrc_init\n");
fprintf(stderr,"     iw = %d\n",iw);
#endif	
		}
 
		/* Option II(b): create Gaussian derivative wavelet to produce effective Ricker
		propagating pulse with specified peak frequency, amplitude at
		spec'd distance. */
		else{
			fprintf(stream,"NOTE: injecting source in the constitutive eq., i.e., source present in press. updates.\n");
			/* RHS wavelet is derivative of gaussian */
			tr->w = igetdgauss(&iw, (m->tsind).dt, fpeak);
			tr->n = 2 * iw+1;
		}
	
		/* source phase - default is zero-phase */
		tr->istart = -iw;
		if (!ps_ffcstring(*par,"waveletphase",&wp)) {
			if (!strcmp(wp,"zerophase")){
				tr->istart = -iw;
			}
			else if (!strcmp(wp,"causal")){
				tr->istart = 0;
			}
/*			else if (!strcmp(wp,"anticausal")) tr->istart=-2*iw;*/
			else {
				fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
				fprintf(stream,"only legit values of waveletphase are \n");
				fprintf(stream,"zerophase and causal\n");
				return E_OTHER;
			}
			fprintf(stream,"NOTE: using wavelet phase: %s\n",wp);
			userfree_(wp);
		}
		else {
			fprintf(stream,"NOTE: wavelet phase not read, using default \"zerophase\"\n");
		}

		/* optionally write out wavelet appearing on RHS */
		tr->idbg = 0;
		ps_ffint(*par, "dump_wavelet", &(tr->idbg));
		
		if (tr->idbg) {
			memcpy(trdbg.data,tr->w,tr->n*sizeof(ireal));
			val.u=1000.0*((m->tsind).dt);
			puthdval(&trdbg,"dt",&val);
			val.h=tr->n;
			puthdval(&trdbg,"ns",&val);
			val.h=((m->tsind).dt)*tr->istart;
			puthdval(&trdbg,"delrt",&val);
		
			if (!(tr->fpdbg=iwave_const_fopen("wavelet.debug","w",NULL,stream))) {
				fprintf(stream,"Error: init_point: failed to open test wavelet file\n");
				return E_FILE;
			}
			else {
				fprintf(stream,"write wavelet trace\n");
			}
			fputtr(tr->fpdbg,&trdbg);
			fflush(tr->fpdbg);
			iwave_fclose(tr->fpdbg);
		}
	}

	
	/* Overall scale factor for source insertion to produce target pulse, per 
	   paper by WWS&TV: contains factors of
	   - 4 pi c^2 dt (from rhs of equation 9, using kappa/rho=c^2 - RHS of 
	   difference scheme is multiplied by kappa * dt);
	   - r (reference distance for normalization, per eqn 13);
	   - reference amplitude;
	   - reciprocal of cell volume, for delta function.
	   Note: for Option Ib (direct insertion of source wavelet, signalled by 
	   refdist<=0) this factor accounts only for time step, bulk modulus, and 
	   cell volume, as this case defines no target scale.
	*/
	prod_d = 1.0;
	for (i = 0; i <  ndim; i++) prod_d *= d[i];
	if (!ps_ffcstring(*par,"source", &srcfile)){
		fprintf(stream,"NOTE: Option I scaling\n");
		if (tr->srcin_flag){
			if (refdist>0) {
				fprintf(stream,"     using scale factor = bou*refdist*refamp*dt/(dxdydz)\n");
				tr->scramp = refbou * refdist * refamp * ((m->tsind).dt) / prod_d;
			}
			else {
				fprintf(stream,"     using scale factor = bou*refamp*dt/(dxdydz)\n");
				tr->scramp = refbou * refamp * ((m->tsind).dt) / prod_d;
			}
		}
		else {
			if (refdist>0) {
				fprintf(stream,"     using scale factor = 4*pi*c^2*refdist*refamp*dt/(dxdydz)\n");
				tr->scramp =  4.0 * 3.1415927 * refvel * refvel * refdist * refamp * ((m->tsind).dt) / prod_d;
			}
			else {
				fprintf(stream,"     using scale factor = kappa*refamp*dt/(dxdydz)\n");
				tr->scramp =  refkappa * refamp * ((m->tsind).dt) / prod_d;
			}
		}
	}
	else {
		fprintf(stream,"NOTE: Option II scaling\n");
		if (tr->srcin_flag){
			if (refdist>0) {
				fprintf(stream,"     using scale factor = bou*refdist*refamp*dt/(dxdydz)\n");
				tr->scramp = refbou * refdist * refamp * ((m->tsind).dt) / prod_d;
			}
			else {
				fprintf(stream,"ERROR: in pointsrc_init, refdist < 0 for case II.\n");
				return 1;
			}
		}
		else {
			if (refdist>0) {
				fprintf(stream,"NOTE: Using scale factor = 4*pi*c^2*refdist*refamp*dt/(dxdydz)\n");
				tr->scramp =  4.0 * 3.1415927 * refvel * refvel * refdist * refamp * ((m->tsind).dt) / prod_d;
			}
			else {
				fprintf(stream,"ERROR: in pointsrc_init, refdist < 0 for case II.\n");
				return 1;
			}
		}

	}

	return 0;  
}

/*----------------------------------------------------------------------------*/
int pointsrc_destroy(POINTSRC * tr) {
/*----------------------------------------------------------------------------*/

#ifdef IWAVE_VERBOSE
	fprintf(stderr,"in pointsrc_destroy\n");
#endif
  	if (tr->w) userfree_(tr->w); 

  	if ( tr->fpsrc ) iwave_fclose(tr->fpsrc);
  	if ( tr->fpdbg ) iwave_fclose(tr->fpdbg);
  	return 0;
}

/*----------------------------------------------------------------------------*/
int pointsrc_run(POINTSRC * tr, IMODEL * m) {
/*----------------------------------------------------------------------------*/

  	int i;

	/* injecting as a body force source, in velocities */
	if (tr->srcin_flag){
		/* key dimn off velocity field - NO-OP if iv!=1 */
//		if ( ((m->tsind).it < tr->n_comp) && ((m->tsind).iv == 1 ) ){
		if ( ((m->tsind).it >= tr->istart) &&
		     ((m->tsind).it <  tr->istart + tr->n) && 
                     ((m->tsind).iv == 1) ){
#ifdef IWAVE_VERBOSE
	fprintf(stderr,">>>> Inside pointsrc_run\n");
	fprintf(stderr,"     evaluating source for iv = %d, it = %d, and istart = %d\n",(m->tsind).iv,(m->tsind).it,tr->istart);
#endif
#ifdef IWAVE_VERBOSE
	fprintf(stderr,"injecting in velocity updates\n");
#endif
    			for (i=0; i < (m->ld_a)._s[D_V0].ndim; i++) {
      				pointsource( tr->is,
		  		             tr->rs,
		  		             tr->order,
		  		             tr->scramp,
		  		             tr->src_d[i] * (tr->w)[ (m->tsind).it - tr->istart ],
       		  	 	             (m->ld_c)._s[D_V[i]],
					     D_V[i] );
    			}
		}
	}
	/* injecting source in pressures */
	else {
		/* key dimn off pressure field - NO-OP if iv!=0 */
// 		if ( ((m->tsind).it < tr->n_comp) && ((m->tsind).iv == 0 ) ){
		if ( ((m->tsind).it >= tr->istart) &&
		     ((m->tsind).it <  tr->istart + tr->n) && 
                     ((m->tsind).iv == 0) ){

    			for (i=0; i < (m->ld_a)._s[D_P0].ndim; i++) {
      				pointsource( tr->is,
		  		             tr->rs,
		  		             tr->order,
		  		             tr->scramp,
		  		             (tr->w)[ (m->tsind).it - tr->istart ],
       		  	 	             (m->ld_c)._s[D_P[i]],
					     D_P[i] );
    			}
		}
  	}
  	return 0;
}

/*----------------------------------------------------------------------------*/
void pointsrc_fprint(POINTSRC const * tr, FILE * fp) {
/*----------------------------------------------------------------------------*/
  	int i;
  	fprintf(fp,"/*---------------------------------------------------------*/\n");
  	fprintf(fp,"POINT SOURCE\n");
	fprintf(fp,"srcin_flag   = %d\n",tr->srcin_flag);
	fprintf(fp,"pulse length = %d\n",tr->n);
  	fprintf(fp,"istart       = %d\n",tr->istart);
  	fprintf(fp,"order        = %d\n",tr->order);
	fprintf(fp,"scramp       = %e\n",tr->scramp);
 	for (i=0;i<RARR_MAX_NDIM;i++)
      		fprintf(fp,"is[%d]=%d rs[%d]=%e src_d[%d]=%e\n", i, tr->is[i], i, tr->rs[i],i,tr->src_d[i]);
    
  	fflush(fp);
}
