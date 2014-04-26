#ifndef __IWAVE_MULTI_SAMPLER__
#define __IWAVE_MULTI_SAMPLER__

#include "sampler.h"

#define FKEY 5

/*MB  Code added by Mario Bencomo for multiple output capabilities.*/
/*****************************************************************************/
/* \section multi_sampler
The MULTI_SAMPLER struct was created for multiple output capabilities. The struct contains an array of samplers along with needed data member for constructing, initializing,running and printing multiple samplers. 
*/
/*****************************************************************************/
typedef struct {
  	char     keyfields[RDOM_MAX_NARR][FKEY]; 	//keys of fields to be sampled
  	SAMPLER  samplers [RDOM_MAX_NARR];		//array of sampler structs
  	int      nsamples; 				//number of fields to be sampled

/* to be assigned by pre-constructor of concrete child class - 
      @param[out] s - (SAMPLER *) sampler to be pre-constructed.
      @param[in] key (const char *) - key for field or field combination, 
      legit values defined in child class. Example: for pressure-velocity 
      acoustics, fields that might be sampled are p and the three velocity 
      components. Accordingly, a child class could assign this pointer to 
      a function which returns appropriate RDOM indices for keys p, v1, v2, v3.

      @return 0 on successful completion, else error code
  */
  	int (*sampler_select)( SAMPLER    * s, 
			       PARARRAY   * pars, 
			       const char * hdrkey, 
			       const char * datakey, 
			       FILE       * stream );
} MULTI_SAMPLER;


/*****************************************************************************/
/* constructor: scans through PARARRAY *pars for data and constructs multiple samples using sampler_select.

    @param[out] m_s - (MULTI_SAMPLER *) multi-sampler object to be constructed.
    @param[in] par - (PARARRAY *) parameter array containing datafile keys etc.
    @param[in] hdrkey - (const char *) key identifying header file name in parameter table 
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_construct( MULTI_SAMPLER * m_s, 
			     PARARRAY      * pars, 
			     const char    * hdrkey, 
			     FILE          * stream  );


/*****************************************************************************/
/* sets members of multi_sampler struct to NULL.

    @param[out] m_s - (MULTI_SAMPLER *) multi-sampler object

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_setnull( MULTI_SAMPLER * m_s );


/*****************************************************************************/
/* initialize multiple samplers.

    @param[out] m_s - (MULTI_SAMPLER *) multi-sampler object to be intialized
    @param[in]  m   - (IMODEL *) model needed for initialization
    @param[in]  par - (PARARRAY *) list of parameters needed for initialization
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_init( MULTI_SAMPLER * m_s, 
                        IMODEL        * m, 
                        PARARRAY      * par, 
                        FILE          * stream );


/*****************************************************************************/
/* destroy multi-sampler object.

    @param[in] m_s - (MULTI_SAMPLER *) multi-sampler object to be destroyed

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_destroy( MULTI_SAMPLER * m_s);


/*****************************************************************************/
/* run multiple samplers

    @param[in] m_s - (MULTI_SAMPLER *) multi-sampler object to be ran
    @param[in] m   - (IMODEL *) model needed for run

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_run( MULTI_SAMPLER * m_s, IMODEL * m);


/*****************************************************************************/
/* print multiple samplers.

    @param[in] m_s - (MULTI_SAMPLER *) multi-sampler object to be printed
    @param[in] stream - (FILE *) verbose output stream

    @return void
 */
/*****************************************************************************/
void multi_sampler_fprint( MULTI_SAMPLER * m_s, FILE * stream);


/*****************************************************************************/
/* set keyfields of multiple samplers. must be done before sampler_select.

    @param[in] m_s - (MULTI_SAMPLER *) multi-sampler object
    @param[in]  par - (PARARRAY *) list of parameters needed for setting keyfields.
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_sampler_set_keyfields( MULTI_SAMPLER * m_s, PARARRAY * pars, FILE * stream );


/*****************************************************************************/
/* writing the traces of multiple samplers

    @param[in] m_s - (MULTI_SAMPLER *) multi-sampler object
    @param[in] d   - (RPNT)
    @param[in] og  - (RPNT)
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
 */
/*****************************************************************************/
int multi_writetraces( MULTI_SAMPLER * m_s, RPNT d, RPNT og, FILE * stream );

#endif