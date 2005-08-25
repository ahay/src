
#ifndef _mt19937ar_h
#define _mt19937ar_h

void init_genrand(unsigned long s);
/* initializes mt[N] with a seed */

double genrand_real1(void);
/* generates a random number on [0,1]-real-interval */

#endif
