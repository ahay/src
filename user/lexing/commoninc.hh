#ifndef _COMMONINC_HH_
#define _COMMONINC_HH_

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


//complex number
#define cpx cpx8

//aux functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }

#define iC(fun)  { int ierr=fun; assert(ierr==0); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<endl; assert(expr); } }

using std::istream;
using std::ostream;
using std::ios_base;
using std::endl;
using std::abs;
using std::cerr;

#endif

