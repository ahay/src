#ifndef __ESG_MULTI_SAMPLER
#define __ESG_MULTI_SAMPLER

#include "multi_sampler.h"
#include "esgn_indices.h"

/** ESG multi-sampler package. */
int esg_sampler_select( SAMPLER * s, 
			PARARRAY * pars, 
			const char * hdrkey, 
			const char * datakey, 
			FILE * stream );

/** ESG pre-constructor */
int esg_multi_sampler_precons( MULTI_SAMPLER * m_s, FILE * stream );

#endif
