#ifndef _sf_bool_h
#define _sf_bool_h

#if defined(__STDC__) && (__STDC_VERSION__ >= 199901L)

#include <stdbool.h>       /* define bool, true, false */

#else

typedef enum {false, true} _Bool;
#define bool _Bool

#endif

#endif
