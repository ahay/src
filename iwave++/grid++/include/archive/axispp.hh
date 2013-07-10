#ifndef __TSOPT_AXIS
#define __TSOPT_AXIS 

#include "contentpackage.hh"
#include <sstream>
#include "parserdecl.hh"
extern "C" {
#include "utils.h"
}

namespace TSOpt {

  using RVL::ContentPackage;
  using RVL::ScalarFieldTraits;

  template<typename Scalar> 
  class Axis {

  public:

    int n;           // number of gridpoints
    Scalar d;        // step
    Scalar o;        // coordinates of origin
    int id;          // axis id number
    Scalar tol;      // comparison tolerance

    Axis()
      : n(1), 
	o(ScalarFieldTraits<Scalar>::Zero()),
	d(ScalarFieldTraits<Scalar>::One()), 
	id(1),
	tol(numeric_limits<Scalar>::epsilon()) {}
    
    bool operator==(Axis<Scalar> const & a) const {
      if (n != a.n) return false;
      if (abs(d-a.d)>tol*abs(d)) return false;
      if (abs(o-a.o)>tol*abs(d)) return false;
      if (id != a.id) return false;
      return true;
    }
    bool operator!=(Axis<Scalar> const & a) const { return !operator==(a); }
  };

}

#endif
