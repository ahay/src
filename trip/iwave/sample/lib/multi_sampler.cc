#include "multi_sampler.h"

#ifndef NSTR
#define NSTR 128
#endif

#define IWAVE_VERBOSE

/*----------------------------------------------------------------------------*/
int multi_sampler_construct( MULTI_SAMPLER * m_s, 
			     PARARRAY      * pars, 
			     const char    * hdrkey, 
			     FILE          * stream  ){
/*----------------------------------------------------------------------------*/
	int err = 0;
	int nsamples = 0;
	char datakey[NSTR];
	char * buff;

	//Counting the number of samples to be outputted, which is computed by sequential scanning of datakeys.
	//NOTE: if datafile keys is missing a sequential key, then the first n sequential keys will be considered.
	while ( !err && nsamples<=RDOM_MAX_NARR ){
		sprintf( datakey, "datafile%d", nsamples+1 );
		err = ps_ffcstring( *pars, datakey, &buff );
		if ( !err && nsamples<RDOM_MAX_NARR ) nsamples++;
#ifdef	IWAVE_VERBOSE
		if (!err) fprintf( stream, "datafile%d was found\n", nsamples );
#endif
	}
	//Case where there were no specified samples found.
	if (nsamples==0){
		fprintf( stream, "Error: Inside multi_sampler_construct\n"
                                 "No type of data output was specified!\n");
		return err;
	}
	m_s->nsamples = nsamples; //setting nsamples

	//Setting datakeys of m_s
	err = multi_sampler_set_keyfields( m_s, pars, stream );
	if (err) {
		fprintf( stream, "Error: in multi_sampler_construct from multi_sampler_setdatakeys\n");
		return err;
	}
	
	//Constructing sampler structs using sampler_select
	int i;
	for (i=0; i<nsamples; i++){
		sprintf( datakey, "datafile%d", i+1 );
		err = m_s->sampler_select( &(m_s->samplers)[i], 
					   pars, 
				     	   hdrkey, 
				     	   datakey, 
				     	   stream );
		if (err){
			fprintf( stream, "Error: in multi_sampler_construct from sampler_select\n"
					 "Invalid value for key %s\n",datakey);
			return err;
		}
	}
	return err;
}

/*----------------------------------------------------------------------------*/
int multi_sampler_setnull( MULTI_SAMPLER * m_s ){
/*----------------------------------------------------------------------------*/
	m_s->nsamples = 0;
	
	int i, j;
	for (i=0; i<RDOM_MAX_NARR; i++){
			snprintf( m_s->keyfields[i], FKEY, "NAN" );
	}
	m_s->sampler_select = NULL;

	return 0;
}

/*----------------------------------------------------------------------------*/
int multi_sampler_init( MULTI_SAMPLER * m_s, 
                        IMODEL        * m, 
                        PARARRAY      * par, 
                        FILE          * stream ){
/*----------------------------------------------------------------------------*/
	int err, i;
	//Initializing sampler structs
	for (i=0; i < m_s->nsamples; i++){
		err = sampler_init( &(m_s->samplers)[i],
		 		    m,
		 		    par, 
		 		    stream);
		if (err) continue;
	}
	if (err) fprintf( stream, "Error: in multi_sampler_init\n");
	return err;

}

/*----------------------------------------------------------------------------*/
int multi_sampler_destroy( MULTI_SAMPLER * m_s){
/*----------------------------------------------------------------------------*/
	int err, i;
	//Destroying sampler structs
	for (i=0; i < m_s->nsamples; i++){
		err = sampler_destroy( &(m_s->samplers)[i]);
		if (err) continue;
	}
	return err;
}

/*----------------------------------------------------------------------------*/
int multi_sampler_run( MULTI_SAMPLER * m_s, IMODEL * m){
/*----------------------------------------------------------------------------*/
	int err, i;
	//Running sampler structs
	for (i=0; i < m_s->nsamples; i++){
		err = sampler_run( &(m_s->samplers)[i], m );
		if (err) continue;
	}
	return err;
}

/*----------------------------------------------------------------------------*/
void multi_sampler_fprint( MULTI_SAMPLER * m_s, FILE * stream) {
/*----------------------------------------------------------------------------*/
	fprintf( stream, "Printing multi-sampler\n" );
	int i;
	for (i=0; i < m_s->nsamples; i++){
		fprintf( stream, "For field %s\n", m_s->keyfields[i] );
		sampler_fprint( &(m_s->samplers)[i], stream );
	}
	return;
}

/*----------------------------------------------------------------------------*/
int multi_sampler_set_keyfields( MULTI_SAMPLER * m_s, PARARRAY * pars, FILE * stream ){
/*----------------------------------------------------------------------------*/
	char datakey[NSTR];
	char * val;
	int i;
	for (i=0; i < m_s->nsamples; i++){
		sprintf( datakey, "datafile%d", i+1 );
		ps_ffcstring( *pars, datakey, &val );
		snprintf( m_s->keyfields[i], FKEY, val );
	}
	return 0;
}
/*----------------------------------------------------------------------------*/
int multi_writetraces( MULTI_SAMPLER * m_s, RPNT d, RPNT og, FILE * stream ){
/*----------------------------------------------------------------------------*/
	int i, err;
	for (i=0; i < m_s->nsamples; i++){
		err=writetraces( &(m_s->samplers[i].t.tg), d, og, stream );
		if (err){
			fprintf(stream,"Error: in multi_writetraces from writetrace\n");
			return err;
		}
	}
	return err;
}