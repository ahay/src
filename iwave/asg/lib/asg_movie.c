#include "asg_movie.h"

int asg_movie_select(const char * key) {
  if (!strcmp(key,"p")) return D_P0;
  if (!strcmp(key,"v0")) return D_V0;
  if (!strcmp(key,"v1")) return D_V1;
  if (!strcmp(key,"v2")) return D_V2;
  return -1;
}

int asg_movie_construct(MOVIE * mt, FILE * stream) {
  movie_setnull(mt);
  mt->iselect=asg_movie_select;
  return 0;
}
