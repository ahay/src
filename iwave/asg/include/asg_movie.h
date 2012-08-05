#ifndef __ASG_MOVIE__
#define __ASG_MOVIE__

#include "movie.h"
#include "sgn_indices.h"

/** ASG movie package. */
int asg_movie_select(const char * key);

/** Movie constructor */
int asg_movie_construct(MOVIE * mt, FILE * stream);

#endif
