/* 
   File: random.c
   --------------
   (Pseudo) random number generator
*/
#include "random.h"

static long seed = 1996;

/* 
   Function: random_init
   ---------------------
   Initialize the seed
*/
void random_init (long iseed)
{
    seed = iseed;
}

/* 
   Function: random0
   -----------------
   Get a random number
   between 0 and 1
*/
float random0 (void)
{
    float rand;
    const int ia = 727, im = 524287;

    seed = (seed*ia)%im;
    rand = ((float) seed - 0.5)/((float) (im - 1));
    return rand;
}
