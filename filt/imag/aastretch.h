#ifndef _aastretch_h
#define _aastretch_h

/* abstract data type */
typedef struct Aamap *aamap;

aamap aastretch_init (int n1, float o1, float d1, int nd);
void aastretch_define (aamap str, 
		       const float *coord, 
		       const float *dt, 
		       const float *amp);
void aastretch_apply (aamap str, const float *ord, float *modl);
void aastretch_close (aamap str);

#endif
