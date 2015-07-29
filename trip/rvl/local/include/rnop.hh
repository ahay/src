/*************************************************************************

Copyright Rice University, 2011.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#ifndef __RVL_RNOP
#define __RVL_RNOP

#include "op.hh"
#include "rnmat.hh"

/** Builds an op (nonlinear function on R^n, or some part of it) out
    of two functions:

    int f(int m, int n, scalar const * in, scalar * out);
    int df(int m, int n, scalar const * in, scalar * jac);

    n is the number of input components, m the number of outputs - so
    f defines a function from R^n to R^m. The dimensions appear as
    arguments to enable sanity-checking, as usual. df has n inputs and m*n
    outputs, namely the jacobian matrix in column major storage.x 
*/

namespace RVL {

  /** FO class encapsulating eval of f */
  template<typename T, int f(int, int, T const *, T *)>
  class CFunction: public BinaryLocalFunctionObject<T> {
  public:
    void operator()(LocalDataContainer<T> & y,
		    LocalDataContainer<T> const & x) {
      int err=0;
      if (err=f(y.getSize(), x.getSize(), x.getData(), y.getData())) {
	RVLException e;
	e<<"Error: CFunction::operator() - err="<<err<<"\n";
	throw e;
      }	
    }
    string getName() const { string tmp = "CFunction"; return tmp; }
  };

  /** FO class encapsulating eval of df. Typical "reaper" design:
      since the thing you want to modify is not a Vector (it's an
      instance of the the matrix LinearOp class GenMat), but it does
      store the data you want to modify in a LocalDataContainer, store
      a reference to this thing on construction. The operator() eval
      method accepts the input LDC as const arg, and does its whatever
      is required to modify your internally stored GenMat. Note that
      the argument(s) of operator() are not modified and are even
      const, so this should be a FOConstEval. You can submit this FO
      to a Vector for evaluation; the Vector will delegate to its DC,
      and if that is an LDC then that will be the argument.
   */ 
  template<typename T, int df(int, int, T const *, T *)>
  class CJacobian: public UnaryLocalFunctionObjectConstEval<T> {
  private:
    GenMat<T> & a;
    CJacobian();
    CJacobian(CJacobian<T, df> const &);
  public:
    CJacobian(GenMat<T> & _a): a(_a) {}
    void operator()(LocalDataContainer<T> const & x) {
      int err=0;
      if (err=df(a.getNRows(), a.getNCols(), x.getData(), a.getData())) {
	RVLException e;
	e<<"Error: CJacobian::operator() - err="<<err<<"\n";
	throw e;
      }
    }
    string getName() const { string tmp = "CJacobian"; return tmp; }
  };

  /** Operator class with matrix Jacobian attribute */
  template<typename T> 
  class OpWithGenMatDeriv: public Operator<T> {
  public:
    virtual GenMat<T> const & getGenMat() const = 0;
  };

  /** Local Operator class based on two C-style functions, one
      computing the function value, the other computing the Jacobian
      matrix.
   */
  template<typename T, int f(int, int, T const *, T *), int df(int, int, T const *, T *)>
  class GenOp: public OpWithGenMatDeriv<T> {
    
  private:

    // domain and range
    RnSpace<T> const & dom;
    RnSpace<T> const & rng;

    // GenMat matrix LinearOp object and intialization flag,
    // modifiable post-construction
    mutable GenMat<T> deriv;
    mutable bool init;

    // common code for intializing the Jacobian (GenMat object),
    // hoisted out. Note that initialization will be done only once,
    // per the design of Operator.
    void buildJacobian(Vector<T> const & x) const {
      try {
	if (!init) {
	  CJacobian<T, df> j(deriv);
	  x.eval(j);
	  init=true;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenOp::buildJacobian\n";
	throw e;
      }
    }

    GenOp();

  protected:
    
    /** \f$y = F(x)\f$ */
    virtual void apply(const Vector<T> & x, 
		       Vector<T> & y) const {
      try {
	CFunction<T, f> ff;
	y.eval(ff,x);
	// why not
	buildJacobian(x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenOp::apply\n";
	throw e;
      }
    }
      
    /** \f$dy = DF(x)dx\f$ */
    virtual void applyDeriv(const Vector<T> & x, 
			    const Vector<T> & dx,
			    Vector<T> & dy) const {
      try {
	buildJacobian(x);
	deriv.applyOp(dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenOp::applyDeriv\n";
	throw e;
      }
    }

    /** \f$dx = DF(x)^*dy\f$ */
    virtual void applyAdjDeriv(const Vector<T> & x, 
			       const Vector<T> & dy,
			       Vector<T> & dx) const {
      try {
	buildJacobian(x);
	deriv.applyAdjOp(dy,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenOp::applyAdjDeriv\n";
	throw e;
      }
    }    

    virtual Operator<T> * clone() const {
      return new GenOp<T, f, df>(dom,rng);
    }

  public:

    /** main constructor: inputs are domain and range RnSpaces */
    GenOp(RnSpace<T> const & _dom, RnSpace<T> const & _rng)
      : dom(_dom), rng(_rng), deriv(dom,rng), init(false) {}
    /** copy constructor - note how important this is, as it's used to
	define the clone method, which is what the OperatorEvaluation
	class uses to create new internal copies of the Operator.
     */
    GenOp(GenOp<T, f, df> const & a)
      : dom(a.dom), rng(a.rng), deriv(dom,rng), init(false) {}
    
    Space<T> const & getDomain() const { return dom; }
    Space<T> const & getRange() const { return rng; }

    /** public access to the internal GenMat - but only if it has been
	initialized */
    GenMat<T> const & getGenMat() const { 
      if (init) return deriv;
      RVLException e;
      e<<"Error: GenOp::getGenMat\n";
      e<<"Jacobian matrix not initialized - must evaluate somewhere first!\n";
      throw e;
    }

    ostream & write(ostream & str) const {
      str<<"GenOp instance\n";
      str<<"domain:\n";
      dom.write(str);
      str<<"range:\n";
      rng.write(str);
      return str; 
    }
  };

}

#endif
