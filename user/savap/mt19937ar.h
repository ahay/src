#ifndef _mt19937ar_h
#define _mt19937ar_h

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

#endif
