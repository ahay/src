#ifndef __TSOPT_SRCWPW__
#define __TSOPT_SRCWPW__

#include "rnspace.hh"
#include "linop.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#else
#include "gridpp.hh"
#include "segypp.hh"
#endif

#include "srcstk.hh"

namespace TSOpt {

  using RVL::RnSpace;
  using RVL::Vector;
  using RVL::LocalVector;
  using RVL::LinearOp;
  using RVL::RVLException;
  using RVL::AssignFilename;
  using RVL::AssignTag;
  using RVL::AssignParams;
  
  /* this class defines a rule for mapping an array of weights to an array 
     source, each trace of which is a reference pulse scaled by a weight. When
     used as source, produces same output as weighted stack of point source
     gathers. */
  
  class SrcWPWOp: public LinearOp<float> {

  private:

    SrcData const & sd;       /* trace data object, creates gather files */
#if IWAVE_USE_MPI
    MPISEGYSpace rng;
#else    
    SEGYSpace rng;            /* space defined by CRG */
#endif
    RnSpace<float> dom;       /* space of weight vectors */
    segy w;                   /* wavelet */

    SrcWPWOp();

  protected:

    LinearOp<float> * clone() const {
      return new SrcWPWOp(*this);
    }

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;
    
  public:
    
    SrcWPWOp(string wavelet,
	     SrcData const & sd);

    SrcWPWOp(SrcWPWOp const & a)
      : sd(a.sd), rng(a.rng), dom(a.dom) {}
    
    ~SrcWPWOp() {}

    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }
    
    ostream & write(ostream & str) const;
  };
}
    
#endif
