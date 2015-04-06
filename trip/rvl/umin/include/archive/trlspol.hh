#ifndef __UMIN_TRLS_POLICY
#define __UMIN_TRLS_POLICY

#include "op.hh"
#include "alg.hh"

namespace RVLUmin {

  using RVL::Vector;
  using RVL::OperatorEvaluation;
  using RVLAlg::Algorithm;

  template<typename T>
  class TRLSSolverPolicy {

    typedef typename ScalarFieldTraits<T>::AbsType AT;

  public:

    virtual Algorithm * build(Vector<T> &,
			      OperatorEvaluation<T> &,
			      AT & rnorm,
			      AT & nrnorm,
			      AT Delta,
			      ostream &str) const = 0;

  };

}
