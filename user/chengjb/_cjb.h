/* cjb.h - include file for general purpose CWP stuff */

#ifndef CJB_H
#define CJB_H

/* INCLUDES */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#ifdef CADDR_T_NOT_DEFINED
typedef char *          caddr_t;
#endif

#ifndef size_t
#define size_t int
#endif
#ifndef NULL
#define NULL	((void *)0)
#endif
#ifndef SEEK_SET
#define SEEK_SET (0)
#endif
#ifndef SEEK_CUR
#define SEEK_CUR (1)
#endif
#ifndef SEEK_END
#define SEEK_END (2)
#endif
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ISIZE sizeof(int)
#define FSIZE sizeof(float)
#define DSIZE sizeof(double)

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef YES
#define YES (1)
#endif
#ifndef NO
#define NO (0)
#endif

#define STREQ(s,t) (strcmp(s,t) == 0)
#define STRLT(s,t) (strcmp(s,t) < 0)
#define STRGT(s,t) (strcmp(s,t) > 0)
#define DIM(a) (sizeof(a)/sizeof(a[0]))
#endif
