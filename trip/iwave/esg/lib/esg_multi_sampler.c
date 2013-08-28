#include "esg_multi_sampler.h"

#ifndef NSTR
#define NSTR 128
#endif

#define IWAVE_VERBOSE
/*----------------------------------------------------------------------------*/
int esg_sampler_select( SAMPLER    * s, 
			PARARRAY   * pars, 
			const char * hdrkey, 
			const char * datakey, 
			FILE       * stream ){
/*----------------------------------------------------------------------------*/
	IPNT sindex;
	RPNT mult, scoord;
	char * val;
	int err = 0;

	//Initializing to 0.
	IASN(sindex,IPNT_0);
	RASN(mult,RPNT_0);
	RASN(scoord,RPNT_0);
	//Extracting key from pars.
	err = ps_ffcstring( *pars, datakey, &val );
	if ( err ){
		fprintf( stream, "Error: in esg_sampler_select:\n"
			         "\"datakey\" = %s was not found!\n",datakey );
		return err;
	}

	if ( !strcmp(val,"p0") ){
		sindex[0] = D_P0;
	          mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_pz_esg.su" );
	}
  	else if ( !strcmp(val,"p1") ){
		sindex[0] = D_P1;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_px_esg.su" );
	}
  	else if ( !strcmp(val,"p2") ){
		sindex[0] = D_P2;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_py_esg.su" );
	}
  	else if ( !strcmp(val,"v0") ){
		sindex[0] = D_V0;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_vz_esg.su" );
	}
  	else if ( !strcmp(val,"v1") ){
		sindex[0] = D_V1;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_vx_esg.su" );
	}
  	else if ( !strcmp(val,"v2") ){
		sindex[0] = D_V2;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_vy_esg.su" );
	}
  	else if ( !strcmp(val,"s0") ){
		sindex[0] = D_S0;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_szx_esg.su" );
	}
  	else if ( !strcmp(val,"s1") ){
		sindex[0] = D_S1;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_sxy_esg.su" );
	}
  	else if ( !strcmp(val,"s2") ){
		sindex[0] = D_S1;
		  mult[0] = 1;
		ps_sfcstring ( *pars, datakey, "trace_szy_esg.su" );
	}
	else{
		fprintf( stream, "Error: %s = %s is not valid for esg!\n", datakey, val );
		return 1;
	}
#ifdef IWAVE_VERBOSE
	ps_ffcstring( *pars, datakey, &val );
	fprintf( stream, "modified in pars: %s = %s\n", datakey, val );
#endif
	//Constructing sampler
	err = sampler_construct( s, pars, sindex, mult, scoord, 0, hdrkey, datakey, stream );
	if ( err ){
		fprintf( stream, "Error: could not construct struct sampler in esg_sampler_select!\n" );
	}
  	return err;
}

/*----------------------------------------------------------------------------*/
int esg_multi_sampler_precons( MULTI_SAMPLER * m_s, PARARRAY * pars, FILE * stream ){
/*----------------------------------------------------------------------------*/
	//Setting multi-sampler m_s to null.
	multi_sampler_setnull( m_s );
	
	//Setting sampler_select of m_s
	m_s->sampler_select = esg_sampler_select;
	return 0;
}