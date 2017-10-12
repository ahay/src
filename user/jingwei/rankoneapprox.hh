#ifndef _RANKONEAPPROX_HH_
#define _RANKONEAPPROX_HH_

#include <rsf.hh>
#include "vecmatop.hh"

using namespace std;
using std::cerr;

int rankoneapprox(const CpxNumMat& M1, const CpxNumMat& M2, CpxNumVec& alpha, CpxNumVec& beta, int npk);

#endif
