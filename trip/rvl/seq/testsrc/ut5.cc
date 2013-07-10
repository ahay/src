#include "seqspace.hh"
#include "multop.hh"
#include "cgnealg.hh"
#include "terminator.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

int main() {

  // cerr<<"coeff list for 1 + x"<<endl;
  list<double> alist;
  alist.push_back(1.0);
  alist.push_back(0.5);
  
  // cerr<<"op of multiplication by 1+x"<<endl;
  PolyMultOp op(alist);

  // cerr<<"vector in domain "<<endl;
  Vector<double> x(op.getDomain());
 
  // cerr<<"vector in range"<<endl;
  Vector<double> b(op.getRange());

  //  cerr<<"assign to b the coeff list of the const poly 1"<<endl;
  list<double> blist;
  blist.push_back(1.0);
  AssignList initb(blist);
  b.eval(initb);

  //  cerr<<"zero out x"<<endl;
  x.zero();

  int nsteps = 50;
  double tol = 1.e-6;
  double rnorm;
  double nrnorm;
  //  cerr<<"create CG step"<<endl;
  CGNEAlg<double> cg(x,op,b,rnorm,nrnorm,tol,tol,nsteps);
  //  cerr<<"run"<<endl; 
  cg.run();

  //  cerr<<"\nComputed Solution:\n";
  x.write(cout);

  //  cerr<<"bye-bye\n";
}
