#ifndef _sf_defs_h
#define _sf_defs_h

#ifndef SF_MAX
#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#endif

#ifndef SF_MIN
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef SF_ABS
#define SF_ABS(a)   ((a) >= 0  ? (a) : (-(a)))
#endif

#ifndef SF_SIG
#define SF_SIG(a)   ((a) >= 0  ?  1  :  -1 )
#endif

#endif
