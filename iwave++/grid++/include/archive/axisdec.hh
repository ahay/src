#ifndef __TSOPT_AXISPP_FUNDEC__
#define __TSOPT_AXISPP_FUNDEC__

/** ContentPackage standalone helper functions for Axis */

#include "axispp.hh"

namespace RVL {

  using TSOpt::Axis;

  template<>
  size_t getDataSize<Axis<float> >(Axis<float> const & ax);

  template<>
  size_t getDataSize<Axis<double> >(Axis<double> const & ax);

  template<>
  ostream & writeMeta<Axis<float> >(Axis<float> const & md, 
				    ostream & s);
  template<>
  ostream & writeMeta<Axis<double> >(Axis<double> const & md, 
				     ostream & s);
}

#endif
