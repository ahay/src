#ifndef __SEQ_FOS
#define __SEQ_FOS

#include "seqdc.hh"

using namespace RVL;

class AssignList: public UnarySeqFO {
private:
  std::list<double> const & inlist;
  AssignList();
public:
  AssignList(std::list<double> const & _inlist): inlist(_inlist) {}
  AssignList(AssignList const & a): inlist(a.inlist) {}
  ~AssignList() {}

  void operator()(SeqDC & x);

  string getName() const { string tmp = "Seq::AssignList"; return tmp; }
};

class ZeroList: public UnarySeqFO {
public:
  ZeroList() {}
  ZeroList(ZeroList const &) {}
  ~ZeroList() {}

  void operator()(SeqDC & x);

  string getName() const { string tmp = "Seq::ZeroList"; return tmp; }
};

class RandList: public UnarySeqFO {
private:
  int minlen;
public:
  RandList(int n=1): minlen(n) {}
  RandList(RandList const & r): minlen(r.minlen) {}
  ~RandList() {}

  void operator()(SeqDC & x);

  string getName() const { string tmp = "Seq::RandList"; return tmp; }
};

class LinCombList: public BinarySeqFO {
private:
  double a, b;
public:
  LinCombList(double _a = 0.0, double _b = 0.0): a(_a), b(_b) {}
  LinCombList(LinCombList const & l): a(l.a), b(l.b) {}
  ~LinCombList() {}

  // y = a*x + b*y
  void operator()(SeqDC & y, SeqDC const & x);
  string getName() const { string tmp = "Seq::LinCombList"; return tmp; }
};

class InnerProdList: public BinarySeqFOR {
public:
  InnerProdList(double bias=0.0): BinarySeqFOR(bias) {}
  InnerProdList(InnerProdList const & ip): BinarySeqFOR(ip.getValue()) {}
  ~InnerProdList() {}
  
  void setValue() { FunctionObjectScalarRedn<double>::setValue(0.0); }
  void setValue(double x) { FunctionObjectScalarRedn<double>::setValue(x); }
  void operator()(SeqDC const & x, SeqDC const & y);
  string getName() const { string tmp = "Seq::InnerProdList"; return tmp; }
};

class PolyMult: public BinarySeqFO {
private:
  SeqDC const & fac;
  PolyMult();
public:
  PolyMult(SeqDC const & _fac): fac(_fac) {}
  PolyMult(PolyMult const & m): fac(m.fac) {}
  ~PolyMult() {}

  // y = fac*x
  void operator()(SeqDC & y, SeqDC const & x);
  string getName() const { string tmp = "Seq::PolyMult"; return tmp; }
};

class PolyMultAdj: public BinarySeqFO {
private:
  SeqDC const & fac;
  PolyMultAdj();
public:
  PolyMultAdj(SeqDC const & _fac): fac(_fac) {}
  PolyMultAdj(PolyMultAdj const & m): fac(m.fac) {}
  ~PolyMultAdj() {}

  // y = fac*x
  void operator()(SeqDC & y, SeqDC const & x);
  string getName() const { string tmp = "Seq::PolyMultAdj"; return tmp; }
};

#endif


