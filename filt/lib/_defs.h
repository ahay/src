#ifndef _sf_defs_h
#define _sf_defs_h

#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))

#define SF_ABS(a)   ((a) >= 0  ? (a) : (-(a)))

#define SF_SIG(a)   ((a) >= 0  ?  1  :  -1 )

#define SF_PI (3.141592653589793)

#endif
