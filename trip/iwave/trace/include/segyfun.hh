#ifndef __CP_SEGY_SPEC__
#define __CP_SEGY_SPEC__

#include "contentpackage.hh"
#include "su.h"
#include "header.h"
#include "segy.h"

/** template specializations for CP wrapper */
namespace RVL {

  template<>
  size_t getDataSize<segy>(segy const & md) {
    size_t n = md.ns;
    return n;
  }

  template<>
  float * newData<float,segy>(segy & md) {
    return md.data;
  }

  template<>
  void deleteData<float,segy>(float ** d,segy ** md) {
    delete *md; *md=NULL; *d=NULL;
  }

  template<>
  ostream & writeMeta<segy>(segy const & md, ostream & s) {
    s<<"SEGY trace\n";
    return s;
  }

}

#endif
