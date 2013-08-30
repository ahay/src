#ifndef __ANSOL_ESG_MOVIE__
#define __ANSOL_ESG_MOVIE__

#include "movie.h"
#include "ansol_esgn.h"

/** ANSOL_ESG movie package. */
int ansol_esg_movie_select(const char * key);

/** Movie constructor */
int ansol_esg_movie_construct(MOVIE * mt, FILE * stream);

#endif