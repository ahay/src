#ifndef __SEQ_POLYMULTOP
#define __SEQ_POLYMULTOP

#include "seqspace.hh"
#include "linop.hh"

class PolyMultOp: public LinearOp<double> {

  /** Defines multiplication by a (fixed) polynomial as a
      RVL::LinearOp. Uses SeqDC wrapper of std::list as coefficient
      sequence data structure, and std::list method to implement. */

private:
  SeqSpace dom;
  SeqDC fac;

  PolyMultOp();

protected:

  LinearOp<double> * clone() const { return new PolyMultOp(*this); }

  void apply(Vector<double> const & x,
	     Vector<double> & y) const;
  void applyAdj(Vector<double> const & x,
		Vector<double> & y) const;

public:

  PolyMultOp(list<double> const & a): dom(), fac() {
    AssignList f(a);
    f(fac);
  }
  PolyMultOp(PolyMultOp const & m): dom(), fac(m.fac) {}
  ~PolyMultOp() {}

  Space<double> const & getDomain() const { return dom; }
  Space<double> const & getRange() const { return dom; }

  ostream & write(ostream & str) const;

};

#endif
  
