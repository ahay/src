#ifndef _COMMONINC_HPP_
#define _COMMONINC_HPP_

//STL stuff
#include <iostream>
#include <fstream>
#include <sstream>

#include <cfloat>
#include <cassert>
#include <cmath>
#include <string>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>

//external libraries
#include "fftw3.h"

#include "blas.h"
#include "lapack.h"

//complex number
typedef std::complex<double> cpx;

//aux functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }

#define iC(fun)  { int ierr=fun; assert(ierr==0); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<endl; assert(expr); } }
//std::cerr<<"wrong"<<std::endl; assert(expr); } }

using std::istream;
using std::ostream;
using std::ios_base;
using std::endl;
using std::abs;
using std::cerr;

#endif

