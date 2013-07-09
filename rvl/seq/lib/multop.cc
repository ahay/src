#include "multop.hh"

void PolyMultOp::apply(Vector<double> const & x,
			 Vector<double> & y) const {
  try {
    PolyMult pm(fac);
    y.eval(pm,x);
  }
  catch (RVLException & e) {
    e<<"\ncalled from PolyMultOp::apply\n";
    throw e;
  }
}

void PolyMultOp::applyAdj(Vector<double> const & x,
			    Vector<double> & y) const {
  try {
    PolyMultAdj pm(fac);
    y.eval(pm,x);
  }
  catch (RVLException & e) {
    e<<"\ncalled from PolyMultOp::applyAdj\n";
    throw e;
  }
}

ostream & PolyMultOp::write(ostream & e) const {
  e<<"PolyMultOp: polynomial multiplication in finite sequence space\n";
  e<<"by factor:\n";
  fac.write(e);
  return e;
}
