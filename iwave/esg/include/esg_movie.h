#ifndef __ESG_MOVIE__
#define __ESG_MOVIE__

#include "movie.h"
#include "esgn_indices.h"

/** ASG movie package. */
int esg_movie_select(const char * key);

/** Movie constructor */
int esg_movie_construct(MOVIE * mt, FILE * stream);

#endif
