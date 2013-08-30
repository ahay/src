#include "ansol_esg_movie.h"

int ansol_esg_movie_select(const char * key) {
	if (!strcmp(key,"ansol_p0")) return D_P0;
  	if (!strcmp(key,"ansol_p1")) return D_P1;
  	if (!strcmp(key,"ansol_p2")) return D_P2;
  	if (!strcmp(key,"ansol_v0")) return D_V0;
  	if (!strcmp(key,"ansol_v1")) return D_V1;
  	if (!strcmp(key,"ansol_v2")) return D_V2;
  	if (!strcmp(key,"ansol_s0")) return D_S0;
  	if (!strcmp(key,"ansol_s1")) return D_S1;
  	if (!strcmp(key,"ansol_s2")) return D_S2;
  	return -1;
}

int ansol_esg_movie_construct(MOVIE * mt, FILE * stream) {
  	movie_setnull(mt);
  	mt->iselect=ansol_esg_movie_select;
  	return 0;
}
