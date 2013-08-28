#ifndef __ANSOL_ESG_MULTI_SAMPLER
#define __ANSOL_ESG_MULTI_SAMPLER

#include "multi_sampler.h"
#include "ansol_esgn.h"

/** ANSOL_ESG multi-sampler package. */
int ansol_esg_sampler_select( SAMPLER * s, 
			      PARARRAY * pars, 
			      const char * hdrkey, 
			      const char * datakey, 
			      FILE * stream );

/** multi-sampler pre-constructor */
int ansol_esg_multi_sampler_precons( MULTI_SAMPLER * m_s, FILE * stream );

#endif
