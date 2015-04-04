#ifndef __ADJSTEPTEST__
#define __ADJSTEPTEST__

#include "acd_gfdm.h"
#include "acd_gfdm2.h"
#include "acd.hh"

using namespace std;

void iwave_rdom_rand(IWAVE *state);
void iwave_rdom_zero(IWAVE *state);
void rdom_rand(RDOM *rd);
void rdom_zero(RDOM *rd);
void iwave_rdom_copy(const IWAVE *state, IWAVE *state1);
void rdom_copy(const RDOM *state, RDOM *state1);
int  ra_a_swap(RARR *arr, RARR *arr1);

// call fw to test if it satisfies adjoint relation
void adjsteptest(std::vector<RDOM *> &iwf, std::vector<RDOM *> &iwa, IWaveInfo const &ic, 
                  void (*tsf)(std::vector<RDOM *> , bool, int, void *fdpars), int,
                  void *fdpars, PARARRAY * pars, FILE * stream);
#endif
