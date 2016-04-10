#ifndef _sf_bool_h
#define _sf_bool_h

#ifndef __cplusplus

#if defined(__STDC__) && (__STDC_VERSION__ >= 199901L)

#include <stdbool.h>       /* define bool, true, false */

#else

typedef enum {false, true} sf_Bool;
#define bool sf_Bool

#endif /* stdbool */

#endif /* c++ */

#endif
