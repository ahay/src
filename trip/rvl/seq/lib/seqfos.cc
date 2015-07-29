#include "seqfos.hh"

using namespace RVL;

void AssignList::operator()(SeqDC & x) {
  std::list<double>::iterator ix = x.get().begin();
  std::list<double>::const_iterator ii = inlist.begin();
  while (ix != x.get().end() && ii != inlist.end()) {
    *ix = *ii;
    ix++;ii++;
  }
  // if input list is shorter than target, pad with zeroes
  while ( ix != x.get().end() ) {
    *ix=0;
    ix++;
  }
  // if input list is longer than target, add remaining 
  // values and lengthen target
  while( ii != inlist.end() ) {
    x.get().push_back(*ii);
    ii++;
  }
}

void ZeroList::operator()(SeqDC & x) {
  std::list<double>::iterator ix = x.get().begin();
  while (ix != x.get().end()) {
    *ix = 0.0;
    ix++;
  }
}

void RandList::operator()(SeqDC & x) {
  std::list<double>::iterator ix = x.get().begin();
  int i = 0;
  while (ix != x.get().end()) {
    *ix = -0.5+ rand()/(RAND_MAX+1.0);
    ix++;i++;
  }
  while (i<minlen) {
    x.get().push_back(-0.5+ rand()/(RAND_MAX+1.0));
    i++;
  }
}

void LinCombList::operator()(SeqDC & y, SeqDC const & x) {
  std::list<double>::iterator iy = y.get().begin();  
  std::list<double>::const_iterator ix = x.get().begin();
  while ((ix != x.get().end()) && (iy != y.get().end())) {
    *iy = a*(*ix) + b*(*iy);
    ix++; iy++;
  }
  // if x summand is longer, treat y as zip after end index of x
  while (ix != x.get().end()) {
    y.get().push_back(a*(*ix));
    ix++;
  }
}

void InnerProdList::operator()(SeqDC const & x, SeqDC const & y) {
  std::list<double>::const_iterator ix = x.get().begin();
  std::list<double>::const_iterator iy = y.get().begin();  
  double val = 0.0;
  while ((ix != x.get().end()) && (iy != y.get().end())) {
    val += (*ix)*(*iy);
    ix++; iy++;
  }
  this->setValue(val);
}

void PolyMult::operator()(SeqDC & y, SeqDC const & x) {

  //  cerr<<"zero output\n";
  ZeroList zip;
  zip(y);
  //  y.write(cerr);

  //  cerr<<"workspace - initialize to x"<<endl;
  SeqDC work;
  AssignList init(x.get());
  init(work);
  //  work.write(cerr);

  //  cerr<<"multiply by factor\n";
  //  fac.write(cerr);

  int deg = 0;
//  std::list<double>::iterator iy = y.get().begin();  
//  std::list<double>::const_iterator ix = x.get().begin();
  std::list<double>::const_iterator ifac = fac.get().begin();
  while (ifac != fac.get().end()) {
    //    cerr<<"degree "<<deg<<endl;
    //    cerr<<"scale shifted workspace by fac coeff = "<<*ifac<<", add to output\n";
    LinCombList lc(*ifac,1.0);
    lc(y,work);

    //    cerr<<"shift workspace to represent multiplication by x^"<<deg<<"\n";
    work.get().push_front(0.0);
    deg++;
    ifac++;
  }
}

void PolyMultAdj::operator()(SeqDC & y, SeqDC const & x) {

  //  cerr<<"start PolyMultAdj\n";
  //  cerr<<"zero output\n";
  ZeroList zip;
  zip(y);
  //  y.write(cerr);

  //  cerr<<"workspace - initialize to x"<<endl;
  SeqDC work;
  AssignList init0(fac.get());
  init0(work);
  AssignList init(x.get());
  init(work);
  //  work.write(cerr);

  //  cerr<<"multiply by factor\n";
  //  fac.write(cerr);

  int deg = 0;
//  std::list<double>::iterator iy = y.get().begin();  
//  std::list<double>::const_iterator ix = x.get().begin();
  std::list<double>::const_iterator ifac = fac.get().begin();
  //  cerr<<"loop\n";
  while (ifac != fac.get().end()) {
    //    cerr<<"degree "<<deg<<endl;
    //    cerr<<"scale shifted workspace by fac coeff = "<<*ifac<<", add to output\n";
    LinCombList lc(*ifac,1.0);
    lc(y,work);

    //    cerr<<"shift workspace to represent division by x^"<<deg<<"\n";
    work.get().pop_front();
    //    cerr<<"loop end\n";
    deg++;
    ifac++;
  }
  //  cerr<<"exit PolyMultAdj\n";
}
