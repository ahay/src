#ifndef _oc_h
#define _oc_h

#include <stdio.h>

#include <rsf.h>

void oc_zero (size_t n, FILE *wall);
void oc_invert(size_t n, FILE *wall);
void oc_dump(size_t n, FILE *wall, sf_file out);

#endif
