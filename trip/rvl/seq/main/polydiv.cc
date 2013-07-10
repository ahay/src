#include "seqspace.hh"
#include "multop.hh"
#include "cgnealg.hh"
#include "terminator.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

int main() {

  cout <<"Polynomial Division in l2 - computing least squares\n";
  cout <<"polynomial approximation to 1/(1+ax) via CG.\n";
  double a;
  cout <<"Enter coeff a of x: ";
  cin  >> a;

  list<double> alist;
  alist.push_back(1.0);
  alist.push_back(a);
  
  // cerr<<"op of multiplication by 1+x"<<endl;
  PolyMultOp op(alist);

  // cerr<<"vector in domain "<<endl;
  Vector<double> x(op.getDomain());
 
  // cerr<<"vector in range"<<endl;
  Vector<double> b(op.getRange());

  // cerr<<"assign to b the coeff list of the const poly 1"<<endl;
  list<double> blist;
  blist.push_back(1.0);
  AssignList initb(blist);
  b.eval(initb);

  // cerr<<"zero out x"<<endl;
  x.zero();

  int nsteps;
  double tol;

  cout<<"Enter limit on number of CG steps: ";
  cin >> nsteps;
  cout<<"Enter tolerance for norm of residual: ";
  cin >> tol;

  if (nsteps < 1) nsteps = 50;
  if (tol < 100*numeric_limits<double>::epsilon()) 
    tol = 100*numeric_limits<double>::epsilon(); 
  // cerr<<"create CG step"<<endl;
  double rnorm;
  double nrnorm;
  CGNEAlg<double> cg(x,op,b,rnorm,nrnorm,tol,tol,nsteps);

  //  cerr<<"run"<<endl; 
  cg.run();

  cout<<"\nComputed Solution:\n";
  x.write(cout);

  // cerr<<"bye-bye\n";
}
