#ifndef _shotfill_h
#define _shotfill_h

void shotfill_init (int nh,float h0, float dh, float ds, float eps);
void shotfill_close (void);
void shotfill_define (float w);
void shotfill_apply (const float complex *s1, const float complex *s2, 
		     float complex *s);

#endif
