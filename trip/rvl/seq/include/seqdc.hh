#ifndef __SEQ_DC
#define __SEQ_DC

#include "contentpackage.hh"
#include "local.hh"
#include "data.hh"

using namespace RVL;

class SeqDC: public DataContainer {

  /** RVL::DataContainer for finite sequences. Wrapper for
      std::list. In addition to usual evaluation methods, provides
      controlled access to internal std::list.

      Evaluation methods admit additional DC arguments only of this
      type. To construct mixed-type evaluations, you will need to
      organize the operations to be performed as a sequence of Unary
      FOs.

      Evaluation is limited to SeqFOs, and to LocalFOs and LocalFORs,
      i.e. function objects whose evaluation methods take
      RVL::LocalDataContainer arguments. In the latter case, the SeqDC
      arguments are buffered to ContentPackages (LDC subtype) and
      passed to the argument LFO or LFOR.

      For definition of SeqFOs, see below.

  */

private:

  std::list<double> datalist;

public:

  // always at least one item
  SeqDC() { datalist.push_back(0.0); }
  SeqDC(SeqDC const & z): datalist(z.datalist) {}
  ~SeqDC() {}

  /** Mutable access to internal std::list data */
  std::list<double> & get() { return datalist; }
  /** Immutable access to internal std::list data */
  std::list<double> const & get() const { return datalist; }

  void eval(FunctionObject & f, 
	    vector<DataContainer const *> & x);
  void eval(FunctionObjectConstEval & f, 
	    vector<DataContainer const *> & x) const;

  ostream & write(ostream & str) const;

};

class UnarySeqFO: public FunctionObject {
  /** This FO type is evaluable only on SeqDCs */
public:
  /** Fundamental evaluation method, used to implement 
      evaluation in SeqDC::eval. */
  virtual void operator()(SeqDC & x) = 0;
};

class BinarySeqFO: public FunctionObject {
  /** This FO type is evaluable only on SeqDCs */
public:
  /** Fundamental evaluation method, used to implement 
      evaluation in SeqDC::eval. */
  virtual void operator()(SeqDC & x, SeqDC const & y) = 0;
};

class UnarySeqFOR: public FunctionObjectScalarRedn<double> {
  /** This FOR type is evaluable only on SeqDCs */
public:
  UnarySeqFOR(double bias=0.0): FunctionObjectScalarRedn<double>(bias) {}
  /** Fundamental evaluation method, used to implement 
      evaluation in SeqDC::eval. */
  virtual void operator()(SeqDC const & x) = 0;
};

class BinarySeqFOR: public FunctionObjectScalarRedn<double> {
  /** This FOR type is evaluable only on SeqDCs */
public:
  BinarySeqFOR(double bias=0.0): FunctionObjectScalarRedn<double>(bias) {}
  /** Fundamental evaluation method, used to implement 
      evaluation in SeqDC::eval. */
  virtual void operator()(SeqDC const & x, 
			  SeqDC const & y) = 0;
};

#endif
