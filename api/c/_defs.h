#ifndef _sf_defs_h
#define _sf_defs_h

#include <float.h>

#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))

#define SF_ABS(a)   ((a) >= 0  ? (a) : (-(a)))
#define SF_SIG(a)   ((a) >= 0  ?  1  :  -1 )

#define SF_NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#define SF_ODD(n)  ((n) & 1)
#define SF_EVEN(n) (!(SF_ODD(n)))

#define SF_PI (3.14159265358979323846264338328)

#define SF_EPS FLT_EPSILON
#define SF_HUGE FLT_MAX

#endif
