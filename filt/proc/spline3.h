#ifndef _spline3_h
#define _spline3_h

void spine3_init(int n1);
void spline3_close (void);
void spline_coeffs(float** table);
float spline_eval(float y);

#endif
