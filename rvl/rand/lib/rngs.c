/* -------------------------------------------------------------------------
 * This is an ANSI C library for multi-stream random number generation.  
 * The use of this library is recommended as a replacement for the ANSI C 
 * rand() and srand() functions, particularly in simulation applications 
 * where the statistical 'goodness' of the random number generator is 
 * important.  The library supplies 256 streams of random numbers; use 
 * SelectStream(s) to switch between streams indexed s = 0,1,...,255.
 *
 * The streams must be initialized.  The recommended way to do this is by
 * using the function PlantSeeds(x) with the value of x used to initialize 
 * the default stream and all other streams initialized automatically with
 * values dependent on the value of x.  The following convention is used 
 * to initialize the default stream:
 *    if x > 0 then x is the state
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively.
 *
 * The generator used in this library is a so-called 'Lehmer random number
 * generator' which returns a pseudo-random number uniformly distributed
 * 0.0 and 1.0.  The period is (m - 1) where m = 2,147,483,647 and the
 * smallest and largest possible values are (1 / m) and 1 - (1 / m)
 * respectively.  For more details see:
 * 
 *       "Random Number Generators: Good Ones Are Hard To Find"
 *                   Steve Park and Keith Miller
 *              Communications of the ACM, October 1988
 *
 * Name            : rngs.c  (Random Number Generation - Multiple Streams)
 * Authors         : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 09-22-98
 * ------------------------------------------------------------------------- 
 */

#include <stdio.h>
#include <time.h>
#include "rngs.h"

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
      
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */


   double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0) 
    seed[stream] = t;
  else 
    seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}


   void PlantSeeds(long x)
/* ---------------------------------------------------------------------
 * Use this function to set the state of all the random number generator 
 * streams by "planting" a sequence of states (seeds), one per stream, 
 * with all states dictated by the state of the default stream. 
 * The sequence of planted states is separated one from the next by 
 * 8,367,782 calls to Random().
 * ---------------------------------------------------------------------
 */
{
  const long Q = MODULUS / A256;
  const long R = MODULUS % A256;
        int  j;
        int  s;

  initialized = 1;
  s = stream;                            /* remember the current stream */
  SelectStream(0);                       /* change to stream 0          */
  PutSeed(x);                            /* set seed[0]                 */
  stream = s;                            /* reset the current stream    */
  for (j = 1; j < STREAMS; j++) {
    x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
    if (x > 0)
      seed[j] = x;
    else
      seed[j] = x + MODULUS;
   }
}


   void PutSeed(long x)
/* ---------------------------------------------------------------
 * Use this function to set the state of the current random number 
 * generator stream according to the following conventions:
 *    if x > 0 then x is the state (unless too large)
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively
 * ---------------------------------------------------------------
 */
{
  char ok = 0;

  if (x > 0)
    x = x % MODULUS;                       /* correct if x is too large  */
  if (x < 0)                                 
    x = ((unsigned long) time((time_t *) NULL)) % MODULUS;              
  if (x == 0)                                
    while (!ok) {
      printf("\nEnter a positive integer seed (9 digits or less) >> ");
      scanf("%ld", &x);
      ok = (0 < x) && (x < MODULUS);
      if (!ok)
        printf("\nInput out of range ... try again\n");
    }
  seed[stream] = x;
}


   void GetSeed(long *x)
/* ---------------------------------------------------------------
 * Use this function to get the state of the current random number 
 * generator stream.                                                   
 * ---------------------------------------------------------------
 */
{
  *x = seed[stream];
}


   void SelectStream(int index)
/* ------------------------------------------------------------------
 * Use this function to set the current random number generator
 * stream -- that stream from which the next random number will come.
 * ------------------------------------------------------------------
 */
{
  stream = ((unsigned int) index) % STREAMS;
  if ((initialized == 0) && (stream != 0))   /* protect against        */
    PlantSeeds(DEFAULT);                     /* un-initialized streams */
}


   void TestRandom(void)
/* ------------------------------------------------------------------
 * Use this (optional) function to test for a correct implementation.
 * ------------------------------------------------------------------    
 */
{
  long   i;
  long   x;
  double u;
  char   ok = 0;  

  SelectStream(0);                  /* select the default stream */
  PutSeed(1);                       /* and set the state to 1    */
  for(i = 0; i < 10000; i++)
    u = Random();
  GetSeed(&x);                      /* get the new state value   */
  ok = (x == CHECK);                /* and check for correctness */

  SelectStream(1);                  /* select stream 1                 */ 
  PlantSeeds(1);                    /* set the state of all streams    */
  GetSeed(&x);                      /* get the state of stream 1       */
  ok = ok && (x == A256);           /* x should be the jump multiplier */    
  if (ok)
    printf("\n The implementation of rngs.c is correct.\n\n");
  else
    printf("\n\a ERROR -- the implementation of rngs.c is not correct.\n\n");
}
