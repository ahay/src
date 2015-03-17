#ifndef __IW_TREE__
#define __IW_TREE__

#include "except.hh"
#include "iwinfo.hh"
#include "iwave.h"

namespace TSOpt {

  using RVL::RVLException;

  // array of IWaveStates with ref object of same
  // type pointing to first half of array
  class IWaveTree {
  private:
    bool own;
    std::vector<IWAVE *> sa;
    std::vector<RDOM *> rd;
    IWaveTree * ref;
    // required for proper destruction
    IWaveInfo const & ic;
    IWaveTree();
    IWaveTree(IWaveTree const &);
    IWaveTree(std::vector<IWAVE *> sv, IWaveInfo const & _ic);

  public:
    IWaveTree(PARARRAY & _pars, FILE * _stream, IWaveInfo const & _ic,
	      int order=0);
    ~IWaveTree();
    
    std::vector<IWAVE *> & getStateArray() { return sa; }
    std::vector<IWAVE *> const & getStateArray() const { return sa; }
    std::vector<IWAVE *> & getRefStateArray() { 
      if (ref) return ref->getStateArray();
      else {
	RVLException e;
	e<<"ERROR: IWaveTree::getRefStateArray()\n";
	e<<"  ref state not initialized - probably order = 0\n";
	e<<"  so no derivative so no reference state\n";
	throw e;
      }
    }
    std::vector<IWAVE *> const & getRefStateArray() const { 
      if (ref) {return ref->getStateArray(); }
      else {
	RVLException e;
	e<<"ERROR: IWaveTree::getRefStateArray()\n";
	e<<"  ref state not initialized - probably order = 0\n";
	e<<"  so no derivative so no reference state\n";
	throw e;
      }
    }
    std::vector<RDOM *> const & getRDOMArray() const { return rd; }
    std::vector<RDOM *> const & getRefRDOMArray() const { return ref->getRDOMArray(); }

    ostream & write(ostream & str) const {
      str<<"IWaveTree, length "<<sa.size()<<"\n";
      return str;
    }
  };
}

#endif
