#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>

#if defined(__sun) || defined(__sun__)
#include <sys/inttypes.h>
#define restrict
#else
#include <stdint.h>
#endif
