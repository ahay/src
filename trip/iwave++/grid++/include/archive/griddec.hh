#ifndef TSOPT_GRIDPP_FUNDEC__
#define TSOPT_GRIDPP_FUNDEC__

#include "gridpp.hh"

namespace RVL {

  template<>
  size_t getDataSize<RARR>(RARR const & g);

  template<>
  ostream & writeMeta<RARR >(RARR const & md, 
			     ostream & s);

  template<>
  ostream & writeMeta<RARR> >(RARR const & md, 
			      ostream & s);

}

#endif
