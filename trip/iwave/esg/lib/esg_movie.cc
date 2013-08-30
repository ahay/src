#include "esg_movie.h"

int esg_movie_select(const char * key) {
  if (!strcmp(key,"p0")) return D_P0;
  if (!strcmp(key,"p1")) return D_P1;
  if (!strcmp(key,"p2")) return D_P2;
  if (!strcmp(key,"v0")) return D_V0;
  if (!strcmp(key,"v1")) return D_V1;
  if (!strcmp(key,"v2")) return D_V2;
  if (!strcmp(key,"s0")) return D_S0;
  if (!strcmp(key,"s1")) return D_S1;
  if (!strcmp(key,"s2")) return D_S2;
  return -1;
}

int esg_movie_construct(MOVIE * mt, FILE * stream) {
  movie_setnull(mt);
  mt->iselect=esg_movie_select;
  return 0;
}
