#ifndef __GRIDPP_FUN__
#define __GRIDPP_FUN__

/** ContentPackage standalone helper functions for Axis */

#include "contentpackage.hh"
#include "rarray.h"

namespace RVL {

  template<>
  size_t getDataSize<RARR>(RARR const & a) {
    size_t n; 
    ra_a_datasize((RARR *)(&a), &n);
    return n;
  }

  template<>
  ireal * newData<ireal,RARR>(RARR & md) {
    int err=0;
    if ((err=ra_allocate(&md))) {
      RVLException e;
      e<<"Error: newData<ireal,RARR> ContentPackage aux fcn\n";
      e<<"from ra_allocate, err="<<err<<"\n";
      throw e;
    }
    return md._s0;
  }

  template<>
  void deleteData<ireal,RARR>(ireal ** d,RARR ** md) {
    //    delete *md; *md=NULL; *d=NULL;
    ra_destroy(*md); delete *md; *md=NULL; *d=NULL;
  }

  template<>
  ostream & writeMeta<RARR>(RARR const & md, ostream & s) {
    s<<"RARR\n";
    return s;
  }

}

#endif

