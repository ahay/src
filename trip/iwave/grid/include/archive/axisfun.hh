#ifndef __AXISPP_FUN__
#define __AXISPP_FUN__

/** ContentPackage standalone helper functions for Axis */

#include "axispp.hh"

namespace RVL {

  using TSOpt::Axis;

  template<>
  size_t getDataSize<Axis<float> >(Axis<float> const & ax) {
    return ax.n;
  }

  template<>
  size_t getDataSize<Axis<double> >(Axis<double> const & ax) {
    return ax.n;
  }

  template<>
  ostream & writeMeta<Axis<float> >(Axis<float> const & md, 
				     ostream & s) {
    s<<"Axis structure, single precision: \n";
    s<<"  number of points = "<<md.n<<"\n";
    s<<"  step             = "<<md.d<<"\n";
    s<<"  origin           = "<<md.o<<"\n";
    s<<"  axis number      = "<<md.id<<"\n";
    return s;
  }

  template<>
  ostream & writeMeta<Axis<double> >(Axis<double> const & md, 
				     ostream & s) {
    s<<"Axis structure, double precision: \n";
    s<<"  number of points = "<<md.n<<"\n";
    s<<"  step             = "<<md.d<<"\n";
    s<<"  origin           = "<<md.o<<"\n";
    s<<"  axis number      = "<<md.id<<"\n";
    return s;
  }

}

#endif
