#ifndef _spline3_h
#define _spline3_h

void spine3_init(int n1);
void spline3_close (void);
void spline_coeffs(float** table);
float spline_eval(float y);

void spine3_init1(int n1, float o1, float d1);
void spline3_close1 (void);
void spline_coeffs1(float* table);
float spline_eval1(float y);

#endif
