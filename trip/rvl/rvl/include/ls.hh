/*************************************************************************

Copyright Rice University, 2004, 2005, 2006
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
#ifndef __LS
#define __LS

#include "op.hh"
#include "functional.hh"

namespace RVL {

  /*template<class Scalar>
  class StdLeastSquaresFcnlGN;
  
  template<class Scalar>
  class RegLeastSquaresFcnlGN;

  template<class Scalar>
  class RegPriorLeastSquaresFcnlGN;
  */

  /** Given an input vector d, this operator implements \f$F(x) = x-d\f$. */
  template<class Scalar>
  class ShiftOperator: public Operator<Scalar> {

  private:

    const Vector<Scalar> & d;
    // disabled
    ShiftOperator();

  protected:

    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	y.copy(x);
	y.linComb(-1.0,d);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ShiftOperator::apply\n";
	throw e;
      }
    }
  
    void applyDeriv(const Vector<Scalar> & x, 
		    const Vector<Scalar> & dx,
		    Vector<Scalar> & dy) const {
      try {
	dy.copy(dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ShiftOperator::applyDeriv\n";
	throw e;
      }
    }

    void applyAdjDeriv(const Vector<Scalar> & x, 
		       const Vector<Scalar> & dy,
		       Vector<Scalar> & dx) const {
      try {
	dx.copy(dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ShiftOperator::applyAdjDeriv\n";
	throw e;
      }
    }
  
    // permits override
    virtual Operator<Scalar> * clone() const {
      return new ShiftOperator<Scalar>(*this);
    }
  
  public:
    /// Usual constructor; just needs a pointer to the linear operator.
    ShiftOperator(const Vector<Scalar> & dd): d(dd) {}
    ShiftOperator(const ShiftOperator<Scalar> & a): d(a.d) {}
    ~ShiftOperator() {}

    /** access to domain, range */
    const Space<Scalar> & getDomain() const { return d.getSpace(); }
    const Space<Scalar> & getRange() const { return d.getSpace(); }

    ostream & write(ostream & str) const {
      str<<"Shift Operator"<<"\n";
      str<<"*** shift vector\n";
      d.write(str);
      return str;
    }
  };

  /** Given an input vector d, and an operator G, this operator
      implements \f$F(x) = G(x)-d\f$. */
  template<class Scalar>
  class ResidualOperator: public Operator<Scalar> {

  private:

    Vector<Scalar> const & d;
    Operator<Scalar> const & G;
    // disabled
    ResidualOperator();

  protected:

    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	this->export_apply(G,x,y);
	y.linComb(-1.0,d);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ResidualOperator::apply\n";
	throw e;
      }
    }
  
    void applyDeriv(const Vector<Scalar> & x, 
		    const Vector<Scalar> & dx,
		    Vector<Scalar> & dy) const {
      try {
	this->export_applyDeriv(G,x,dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ResidualOperator::applyDeriv\n";
	throw e;
      }
    }

    void applyAdjDeriv(const Vector<Scalar> & x, 
		       const Vector<Scalar> & dy,
		       Vector<Scalar> & dx) const {
      try {
	this->export_applyAdjDeriv(G,x,dy,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ResidualOperator::applyAdjDeriv\n";
	throw e;
      }
    }
  
    // permits override
    virtual Operator<Scalar> * clone() const {
      return new ResidualOperator<Scalar>(*this);
    }
  
  public:
    /// Usual constructor; just needs a pointer to the linear operator.
    ResidualOperator(Operator<Scalar> const & GG,
		     Vector<Scalar> const & dd): d(dd), G(GG) {
      if (G.getRange() != d.getSpace()) {
	RVLException e;
	e<<"Error: ResidualOperator constructor\n";
	e<<"input vector not in range of operator\n";
	e<<"input vector:\n";
	d.write(e);
	e<<"operator:\n";
	G.write(e);
	throw e;
      }
    }
      ResidualOperator(ResidualOperator<Scalar> const & a): d(a.d), G(a.G) {}
    ~ResidualOperator() {}

    /** access to domain, range */
    const Space<Scalar> & getDomain() const { return G.getDomain(); }
    const Space<Scalar> & getRange() const { return G.getRange(); }

    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      return G.getMaxStep(x,dx);
    }

    ostream & write(ostream & str) const {
      str<<"Residual Operator"<<"\n";
      str<<"*** data vector:\n";
      d.write(str);
      str<<"*** operator:\n";
      G.write(str);
      return str;
    }
  };

  /** This functional is the standard Euclidean Form \f$f(x) = 0.5
      |x|^{2}\f$. Note that this form is differentiable over the
      reals, and the gradient provided is its derivative. However it
      is actually not differentiable for complex Scalar types, and
      accordingly the derivative test would fail.

      A call to testRealOnly is included to catch this point at
      compile time.
  */
  template<class Scalar>
  class EuclideanForm: public Functional<Scalar> {

  private:

    const Space<Scalar> & sp;
    EuclideanForm();

  protected: 

    void apply(const Vector<Scalar> & x,
	       Scalar & val) const {
      try {
	if (sp != x.getSpace()) {
	  RVLException e;
	  e<<"Error: EuclideanForm::apply\n";
	  e<<"input vector not in domain\n";
	  throw e;
	}
	val = 0.5*x.inner(x);
      }
      catch(RVLException & e) {
	e<<"\ncalled from EuclideanForm::apply\n";
	throw e;
      }
    }

    void applyGradient(const Vector<Scalar> & x,
		       Vector<Scalar> & g) const {
      try {
	if (sp != x.getSpace()) {
	  RVLException e;
	  e<<"Error: EuclideanForm::applyGradient\n";
	  e<<"input vector not in domain\n";
	  throw e;
	}
	g.copy(x);
      }
      catch(RVLException & e) {
	e<<"\ncalled from EuclideanForm::applyGradient\n";
	throw e;
      }
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & delx,
		      Vector<Scalar> & dely) const {
      try {
	if (sp != delx.getSpace()) {
	  RVLException e;
	  e<<"Error: EuclideanForm::applyHessian\n";
	  e<<"input vector not in domain\n";
	  throw e;
	}
	dely.copy(delx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from EuclideanForm::applyGradient\n";
	throw e;
      }
    }

    // permits override
    virtual Functional<Scalar> * clone() const { 
      return new EuclideanForm<Scalar>(*this); 
    }

  public:

    EuclideanForm(const Space<Scalar> & _sp): sp(_sp) { testRealOnly<Scalar>(); }
    EuclideanForm(const EuclideanForm<Scalar> & q): sp(q.sp) { testRealOnly<Scalar>(); }
    ~EuclideanForm() {}

    // access to domain
    const Space<Scalar> & getDomain() const { return sp; }

    ostream & write(ostream & str) const {
      str<<"Euclidean quadratic form - length squared\n";
      return str;
    }
  };

  /** QuadraticForm creates a function of the form
      \f$
      x \mapsto \frac{1}{2} |A x|^2
      \f$
      in which A is a linear operator.

      Only differentiable when Scalar is not a complex type.

  */
  template<class Scalar>
  class QuadraticForm: public Functional<Scalar> {

  private:

    const LinearOp<Scalar> & A;
    mutable Vector<Scalar> r;
    mutable bool applied;

    QuadraticForm();

  protected: 

    void apply(const Vector<Scalar> & x,
	       Scalar & val) const {
      try {
	if (!applied) {
	  A.applyOp(x,r);
	  applied = true;
	}
	val =  0.5*r.normsq();;
      }
      catch(RVLException & e) {
	e<<"\ncalled from QuadraticForm::apply\n";
	throw e;
      }
    }

    void applyGradient(const Vector<Scalar> & x,
		       Vector<Scalar> & g) const {
      try {
	if (!applied) {
	  A.applyOp(x,r);
	  applied = true;
	}
	A.applyAdjOp(r,g);
      }
      catch(RVLException & e) {
	e<<"\ncalled from QuadraticForm::applyGradient\n";
	throw e;
      }
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & delx,
		      Vector<Scalar> & dely) const {
      try {
	Vector<Scalar> z(A.getRange());
	A.applyOp(delx,z);
	A.applyAdjOp(z,dely);
      }
      catch (RVLException & e) {
	e<<"\ncalled from QuadraticForm::applyGradient\n";
	throw e;
      }
    }

    // permits override
    virtual Functional<Scalar> * clone() const { 
      return new QuadraticForm<Scalar>(*this); 
    }

  public:

    QuadraticForm(const LinearOp<Scalar> & AA)
      : A(AA), r(A.getDomain()), applied(false) { testRealOnly<Scalar>(); }

    QuadraticForm(const QuadraticForm<Scalar> & q)
      : A(q.A), r(A.getDomain()), applied(false) { testRealOnly<Scalar>(); }
		
    ~QuadraticForm() {}

    const Space<Scalar> & getDomain() const { return A.getDomain(); }

    ostream & write(ostream & str) const {
      str<<"Quadratic form with operator\n";
      A.write(str);
      return str;
    }
  };

  /** ShiftedQuadraticForm creates a function of the form
      \f$
      x \mapsto \frac{1}{2} |A x - b|^2
      \f$
      in which A is a linear operator, b a vector.

      Only differentiable when Scalar is not a complex type.

  */
  template<class Scalar>
  class ShiftedQuadraticForm: public Functional<Scalar> {

  private:

    const LinearOp<Scalar> & A;
    const Vector<Scalar> & b;
    mutable Vector<Scalar> r;
    mutable bool applied;

    ShiftedQuadraticForm();

  protected: 

    void apply(const Vector<Scalar> & x,
	       Scalar & val) const {
      try {
	if (!applied) {
	  A.applyOp(x,r);
	  r.linComb(-1.0,b);
	  applied = true;
	}
	val = 0.5*r.normsq();
      }
      catch(RVLException & e) {
	e<<"\ncalled from QuadraticForm::apply\n";
	throw e;
      }
    }

    void applyGradient(const Vector<Scalar> & x,
		       Vector<Scalar> & g) const {
      try {
	if (!applied) {
	  A.applyOp(x,r);
	  r.linComb(-1.0,b);
	  applied = true;
	}
	A.applyAdjOp(r,g);
      }
      catch(RVLException & e) {
	e<<"\ncalled from QuadraticForm::applyGradient\n";
	throw e;
      }
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & delx,
		      Vector<Scalar> & dely) const {
      try {
	Vector<Scalar> z(A.getRange());
	A.applyOp(delx,z);
	A.applyAdjOp(z,dely);
      }
      catch (RVLException & e) {
	e<<"\ncalled from QuadraticForm::applyGradient\n";
	throw e;
      }
    }

    // permits override
    virtual Functional<Scalar> * clone() const { 
      return new ShiftedQuadraticForm<Scalar>(*this); 
    }

  public:

    ShiftedQuadraticForm(const LinearOp<Scalar> & AA,
			 const Vector<Scalar> & bb)
      : A(AA), b(bb), r(A.getRange()), applied(false) { testRealOnly<Scalar>(); }

    ShiftedQuadraticForm(const ShiftedQuadraticForm<Scalar> & q)
      : A(q.A), b(q.b), r(A.getRange()), applied(false) { testRealOnly<Scalar>(); }
		
    ~ShiftedQuadraticForm() { }

    const Space<Scalar> & getDomain() const { return A.getDomain(); }

    ostream & write(ostream & str) const {
      str<<"ShiftedQuadratic form with linear operator\n";
      A.write(str);
      str<<"and shift vector\n";
      b.write(str);
      return str;
    }
  };

  /** LeastSquaresFcnlGN creates a least squares objective
      function from an operator. This is a so called \em bridge
      class, i.e. it exists merely to make objects of some class
      behave like objects of another class. The input object is an
      operator \f$ A \f$ (an instance of Operator ). The class
      uses a reference to this object to construct the least squares
      functinal \f$ x \mapsto \frac{1}{2} \|A[x]\|^2 \equiv
      J[x] \f$, its gradient \f$ x \mapsto DA[x]^T(A[x]) =
      \nabla J[x] \f$, and the Gauss-Newton approximation to its
      Hessian operator \f$ x \mapsto DA[x]^TDA[x] \f$ and
      presents the result as an instance of the Functional
      interface.

      As with all LS functions, only defines a differentiable function
      is Scalar is a real type.
  */

  template<class Scalar>
  class LeastSquaresFcnlGN: public Functional<Scalar> {

  private:

    EuclideanForm<Scalar> sql;
    FcnlOpComp<Scalar> work;

    // default constructor--disabled
    LeastSquaresFcnlGN();

  protected:
    
    void apply(Vector<Scalar> const & x,
	       Scalar & val) const {
      try {
	this->export_apply(work,x,val);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LeastSquaresFcnlGN::apply\n";
	throw e;
      }
    }
    void applyGradient(Vector<Scalar> const & x,
		       Vector<Scalar> & g) const {
      try {
	this->export_applyGradient(work,x,g);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LeastSquaresFcnlGN::applyGradient\n";
	throw e;
      }
    }
    void applyHessian(Vector<Scalar> const & x,
		      Vector<Scalar> const & dx,
		      Vector<Scalar> & dy) const {
      try {
	this->export_applyHessian(work,x,dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LeastSquaresFcnlGN::applyHessian\n";
	throw e;
      }
    }

    virtual Functional<Scalar> * clone() const { 
      return new LeastSquaresFcnlGN<Scalar>(*this); 
    }

  public:  
  
    /** Usual constructor */
    LeastSquaresFcnlGN(Operator<Scalar> const & op)
      : sql(op.getRange()), work(sql,op) {
      testRealOnly<Scalar>(); 
    }
    /** Copy constructor */
    LeastSquaresFcnlGN(const LeastSquaresFcnlGN<Scalar> & J)
      : sql(J.sql), work(J.work) { testRealOnly<Scalar>(); }
  
    // Destructor.
    virtual ~LeastSquaresFcnlGN() {}

    Space<Scalar> const & getDomain() const { return work.getDomain(); }
 
    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      try {
	return work.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LeastSquaresFcnlGN::getMaxStep\n";
	throw e;
      }
    }
    
    ostream & write(ostream & str) const {
      str<<"Least-Squares Gauss-Newton functional: expressed as = \n";
      work.write(str);
      return str;
    }
  };

  /** StdLeastSquaresFcnlGN creates a least squares objective function
      from an operator and a data vector. This is a so called \em
      bridge class, i.e. it exists merely to make objects of some
      class behave like objects of another class. The input objects
      are an operator \f$ A \f$ (an instance of Operator) and a vector
      \f$ d \f$ (an instance of Vector). The class uses a reference to
      this object to construct the least squares functinal \f$ x
      \mapsto \frac{1}{2} \|A[x]-d\|^2 \equiv J[x] \f$, its gradient
      \f$ x \mapsto DA[x]^T(A[x]-d) = \nabla J[x] \f$ and the
      Gauss-Newton approximation to its Hessian operator \f$ x \mapsto
      DA[x]^TDA[x] \f$ and presents the result as an instance of the
      Functional interface.

      Defines differentiable function only if Scalar is a real type, 
      so this property is tested at compile time.
  */

  template<class Scalar>
  class StdLeastSquaresFcnlGN: public Functional<Scalar> {

  private:

    // Internal variables
    EuclideanForm<Scalar> sql;
    ResidualOperator<Scalar> res;
    FcnlOpComp<Scalar> work;

    // default constructor--disabled
    StdLeastSquaresFcnlGN();

  protected:
    
    void apply(Vector<Scalar> const & x,
	       Scalar & val) const {
      try {
	this->export_apply(work,x,val);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdLeastSquaresFcnlGN::apply\n";
	throw e;
      }
    }
    void applyGradient(Vector<Scalar> const & x,
		       Vector<Scalar> & g) const {
      try {
	this->export_applyGradient(work,x,g);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdLeastSquaresFcnlGN::applyGradient\n";
	throw e;
      }
    }
    void applyHessian(Vector<Scalar> const & x,
		      Vector<Scalar> const & dx,
		      Vector<Scalar> & dy) const {
      try {
	this->export_applyHessian(work,x,dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdLeastSquaresFcnlGN::applyHessian\n";
	throw e;
      }
    }

    virtual Functional<Scalar> * clone() const { 
      return new StdLeastSquaresFcnlGN<Scalar>(*this); 
    }
  
  public:  
  
    /** Usual constructor */
    StdLeastSquaresFcnlGN(Operator<Scalar> const & oper, 
			  Vector<Scalar> const & d)
      : sql(oper.getRange()), 
	res(oper,d),
        work(sql,res) {
      testRealOnly<Scalar>(); 
    }

    /** Copy constructor */
    StdLeastSquaresFcnlGN(const StdLeastSquaresFcnlGN<Scalar> & JJ)
      : sql(JJ.sql), 
	res(JJ.res), 
	work(JJ.work) {
      testRealOnly<Scalar>();
    }
  
    // Destructor.
    virtual ~StdLeastSquaresFcnlGN() {}

    Space<Scalar> const & getDomain() const { return work.getDomain(); }
 
    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      try {
	return res.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdLeastSquaresFcnlGN::getMaxStep\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"Standard Least-Squares Gauss-Newton functional\n";
      str<<"based on ResidualOperator\n";
      res.write(str);
      return str;
    }
  };

  /** RegLeastSquaresFcnlGN creates a least squares objective function
      with a penalty term from an operator, a linear operator, a data
      vector, and a prior. This is a so called \em bridge class,
      i.e. it exists merely to make objects of some class behave like
      objects of another class. The input objects are an operator A
      (an instance of Operator), a linear operator R (an instance of
      LinearOp), and vectors d and r0 (instances of Vector). Note that
      r0 is in the range of R - if the application provides a model
      space prior x0, then r0=Rx0 must be precomputed and passed to
      the constructor. The class uses references to these objects to
      construct the least squares functional \$f x \mapsto \frac{1}{2}
      (\|A[x]-d\|^2 + \lambda^2 \|Rx-r_0)\|^2) \equiv J[x] \f$, its
      gradient \f$ x \mapsto DA[x]^T(A[x]-d) + \lambda^2 R^T(Rx-r_0) =
      \nabla J[x] \f$, and the Gauss-Newton approximation to its
      Hessian operator \f$ x \mapsto DA[x]^TDA[x] + \lambda^2 R^TR
      \f$, and presents the result as an instance of the Functional
      interface.

      Real scalar types only.
  */

  template<class Scalar>
  class RegLeastSquaresFcnlGN: public Functional<Scalar> {

  private:

    // default constructor--disabled
    RegLeastSquaresFcnlGN();
    StdLeastSquaresFcnlGN<Scalar> f1;
    ShiftedQuadraticForm<Scalar> f2;
    Scalar lambda;
    LinCombFunctional<Scalar> work;

  protected:

    void apply(Vector<Scalar> const & x,
	       Scalar & val) const {
      try {
	this->export_apply(work,x,val);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RegLeastSquaresFcnlGN::apply\n";
	throw e;
      }
    }
    void applyGradient(Vector<Scalar> const & x,
		       Vector<Scalar> & g) const {
      try {
	this->export_applyGradient(work,x,g);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RegLeastSquaresFcnlGN::applyGradient\n";
	throw e;
      }
    }
    void applyHessian(Vector<Scalar> const & x,
		      Vector<Scalar> const & dx,
		      Vector<Scalar> & dy) const {
      try {
	this->export_applyHessian(work,x,dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RegLeastSquaresFcnlGN::applyHessian\n";
	throw e;
      }
    }

    virtual Functional<Scalar> * clone() const { 
      return new RegLeastSquaresFcnlGN<Scalar>(*this); 
    }

  public:  
  
    /** Usual constructor */
    RegLeastSquaresFcnlGN(Operator<Scalar> const & oper, 
			  LinearOp<Scalar> const & reg,
			  Vector<Scalar> const & d,
			  Vector<Scalar> const & r0,
			  Scalar lam)
      : f1(oper,d), f2(reg,r0), lambda(lam), work(1.0,f1,lam,f2) {
      try {
	testRealOnly<Scalar>(); 
      }
      catch (RVLException & e) {
	e<<"\ncalled from RegLeastSquaresFcnlGN constructor\n";
	throw e;
      }
    }

    /** Copy constructor */
    RegLeastSquaresFcnlGN(const RegLeastSquaresFcnlGN<Scalar> & rls)
      : f1(rls.f1), f2(rls.f2), lambda(rls.lambda), work(rls.work) { 
      testRealOnly<Scalar>(); 
    }
  
    // Destructor.
    virtual ~RegLeastSquaresFcnlGN() {}

    Space<Scalar> const & getDomain() const { return work.getDomain(); }

    ostream & write(ostream & str) const {
      str<<"Regularized Least-Squares Gauss-Newton functional. realized as\n";
      work.write(str);
      return str;
    }
  };

    
    
    
  /** Implements bounds test by returning infty if bounds are
      violated.  The data member b should return true if x is
      feasible, else false - it is an instance of RVL::Oracle. Note
      that feasibility need not mean violation of bounds!
  */
  
  template<typename Scalar>
  class FunctionalBd: public Functional<Scalar> {

  private:

    //    Functional<Scalar> const & f;
    Functional<Scalar> * f;
    RVL::Oracle<Vector<Scalar> > const & b;

  protected:

    void apply(const Vector<Scalar> & x, 
	       Scalar & val) const {
      try {
	val=numeric_limits<Scalar>::max();
	if (!b.isFeasible(x)) return;
	RVL::Functional<Scalar>::export_apply(*f,x,val);
      }
      catch (RVLException & e) {
	val=numeric_limits<Scalar>::max();
      }
    }

    void applyGradient(const Vector<Scalar> & x, 
		       Vector<Scalar> & g) const {
      try {
	if (!b.isFeasible(x)) {
	  RVLException e;
	  e<<"Error: FuncionalBd::applyGradient\n";
	  e<<"infeasible point\n";
	  throw e;
	}
	RVL::Functional<Scalar>::export_applyGradient(*f,x,g);
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalBd::applyGradient\n";
	throw e;
      }
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx, 
		      Vector<Scalar> & dy) const {
      try {
	if (!b.isFeasible(x)) {
	  RVLException e;
	  e<<"Error: FuncionalBd::applyHessian\n";
	  e<<"infeasible point\n";
	  throw e;
	}
	RVL::Functional<Scalar>::export_applyHessian(*f,x,dx,dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalBd::applyHessian\n";
	throw e;
      }
      
    }
    
    // permits override
    virtual Functional<Scalar> * clone() const { 
      return new FunctionalBd(*this); 
    }
    
  public:
    
    FunctionalBd(Functional<Scalar> const & _f,
		 RVL::Oracle<Vector<Scalar> > const & _b) 
      : f(NULL), b(_b) { testRealOnly<Scalar>(); this->export_clone(_f,&f); }

    FunctionalBd(FunctionalBd<Scalar> const & a)
      :	f(NULL), b(a.b) { testRealOnly<Scalar>(); this->export_clone(*(a.f),&f); }

    ~FunctionalBd() { if (f) delete f; }

    Space<Scalar> const & getDomain() const { return f->getDomain(); }
    Scalar getMaxStep(Vector<Scalar> const & x, 
		      Vector<Scalar> const & dx) const { 
      return f->getMaxStep(x,dx);
    }

    Functional<Scalar> const & getFunctional() const { return *f; }

    ostream & write(ostream & str) const {
      str<<"FunctionalBd: functional returning infty at infeasible points\n";
      str<<"data member Functional:\n";
      f->write(str);
      str<<"bounds test:\n";
      b.write(str);
      return str;
    }

  };

  /** Convenient definition of typical upper and lower bound constraints,
      often (though not necessarily) implemented componentwise. */
  
  template<typename Scalar>
  class ULBoundsTest: public RVL::Oracle<Vector<Scalar> > {
    
  private: 
    Vector<Scalar> const & ub;
    Vector<Scalar> const & lb;
    FunctionObjectScalarRedn<Scalar> & minfo;

  public:

    ULBoundsTest(Vector<Scalar> const & _lb,
		 Vector<Scalar> const & _ub,
		 FunctionObjectScalarRedn<Scalar> & _minfo)
      : ub(_ub), lb(_lb), minfo(_minfo) {
      testRealOnly<Scalar>(); 
      if (lb.getSpace() != ub.getSpace()) {
	RVLException e;
	e<<"Error: ULBoundsTest constructor\n";
	e<<"upper, lower bound Vectors not in same Space\n";
	e<<"lower bd vector:\n";
	lb.write(e);
	e<<"upper bd vector:\n";
	ub.write(e);
	throw e;
      }
      minfo.setValue();
    }

    ULBoundsTest(ULBoundsTest<Scalar> const & t)
      : lb(t.lb), ub(t.ub), minfo(t.minfo) {
      testRealOnly<Scalar>(); 
      minfo.setValue();
    }		

    ~ULBoundsTest() {}

    bool isFeasible(Vector<Scalar> const & x) const {
      // this function should be called only once per x, so
      // allocate workspace locally
      try {
	if (x.getSpace() != lb.getSpace()) return false; 
	Vector<Scalar> work(x.getSpace());
	Components<Scalar> cwork(work);
	// test for bound violation, return false if found
	float one = ScalarFieldTraits<float>::One();
	float zip = ScalarFieldTraits<float>::Zero();
	int size = cwork.getSize();
	float valmin = zip;
	minfo.setValue(zip);
	work.copy(x);
	work.linComb(one,ub,-one);
	for(int i=0; i< size; i++) {
	  cwork[i].eval(minfo);
	  // cerr<<"min(ub - x)["<<i<<"] ="<<minfo.getValue()<<endl;
	  if (minfo.getValue() < valmin)
	    valmin = minfo.getValue();
	}
	// cerr<<"min(ub - x) = valmin ="<<valmin<<endl;
	work.copy(x);
	work.linComb(-one,lb);
	for(int i=0; i<size; i++){
	  cwork[i].eval(minfo);
	  // cerr<<"min(ub - x, x - lb)["<<i<<"]="<<minfo.getValue()<<endl;
	  if (minfo.getValue() < valmin)
	    valmin = minfo.getValue();
	}
	// cerr<<"min(ub - x, x - lb) = valmin ="<<valmin<<endl;
	if (valmin<zip) return false;
	// passed all tests
	return true;
      }
      catch (RVLException & e) {
	return false;
      }
    }

    ostream & write(ostream & str) const {
      str<<"Standard upper/lower bounds test\n";
      str<<"Lower bound vector:\n";
      lb.write(str);
      str<<"upper bound vector:\n";
      ub.write(str);
      str<<"min function:\n";
      minfo.write(str);
      return str;
    }
  };
}
#endif
