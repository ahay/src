#ifndef _randqm_h
#define _random_h

/* 
   File: random.h
   --------------
   (Pseudo) random number generator
*/

/* 
   Function: random_init
   ---------------------
   Initialize the seed
*/
void random_init (long iseed);

/* 
   Function: random0
   -----------------
   Get a random number
   between 0 and 1
*/
float random0 (void);

#endif

/* 	$Id: random.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
