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

#ifndef __RVL_LINOP
#define __RVL_LINOP

#include "space.hh"
#include "write.hh"

namespace RVL {

  // because of this forward declaration, this header file MUST be
  // referenced only through inclusion in op.hh
  template<typename Scalar>
  class Operator;

  /** Principal RVL interface for Linear Operators. Since linear
      operators in Hilbert space have well-defined adjoints, this
      class really defines a pair of operators, adjoint to each
      other.

      The abstract public interface consists of 

      <ul>
      
      <li>getDomain, getRange: statements of domain and range</li>
      
      <li>applyOp, applyAdjOp: application of the operator and its adjoint</li>
      
      <ul>

      The getDomain() method returns a const reference to the
      RVL::Space representing the domain; getRange() does likewise for
      the range. These are pure virtual methods, so will need to be
      implemented in instantiable derived classes. Typical
      implementations of getDomain and getRange will expose references
      to const Space data members.

      The applyOp() method accepts a const RVL::Vector argument
      representing the input vector, and a mutable RVL::Vector
      argument representing the ouput. The input must be a member of
      the domain, and the output a member of the range, and these
      conditions are checked in the method body using
      RVL::Space::operator=, RVL::Vector::getSpace, and the getDomain
      and getRange methods of this class. Similarly the applyAdjOp
      method implements application of the adjoint, and sanity-checks
      its arguments as well. Note that both applyOp and applyAdjOp are
      implemented, and may not be overridden in child classes.

      The applyOp and applyAdjOp methods delegate the computations
      inherent in the action of the linear operator and its adjoint to
      two pure virtual protected methods, apply and applyAdj, which
      accept the same arguments. The subclass writer will thus need to
      implement apply and applyAdj. This device permits the sanity
      tests for membership in domain and range to be hard-coded in
      applyOp and applyAdjOp, relieving the subclass writer of the
      necessity to test the dimensional and other conditions implicit
      in membership.

      A protected pure virtual copy constructor (clone) must also be
      supplied; a typical implementation delegates to the copy
      constructor of the subclass.

      In summary, the subclass writer must implement three protected methods:
      <ul>
      <li>clone: usually calls copy constructor</li>

      <li>apply: code for action of operator, can assume membership of
      input and output in domain and range without further tests</li>

      <li>applyAdj: code for adjoint action of operator, can assume
      membership of input and output in domain and range without
      further tests</li>

      </ul>

      and two public methods (usually by returning reference member data):

      <ul>

      <li>getDomain: returns const reference to domain RVL::Space</li>

      <li>getRange: returns const reference to range RVL::Space</li>

      </ul>

      A unit test for validity of the adjoint relationship is supplied
      as part of RVL, in the form of a standalone function
      (AdjointTest). It is HIGHLY RECOMMENDED that every concrete
      LinearOp implementation be subjected to this test.
  */

  template <class Scalar>
  class LinearOp: public Operator<Scalar> {
    
  protected:

    /** operator new - not available to general public, but
	available to children who will use it to define a 
	clone method, for example 

	Version 1.0: user control, for the bold. "There are old
	sailors, and there are bold sailors, but there are no old,
	bold sailors."
    */

#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

    // note: clone and apply are inherited from Operator

    /** Evaluation of adjoint linear operator on constant input vector
	x, output written on mutable vector y. Accessed only through
	public applyOp method, so creator of subclasses may assume
	that input is in domain, output is in range */
    virtual void applyAdj(const Vector<Scalar> & x,
			  Vector<Scalar> & y) const = 0;

    /** implemented derivative method - a linear operator is its own derivative */
    void applyDeriv(const Vector<Scalar> & x, 
		    const Vector<Scalar> & dx,
		    Vector<Scalar> & dy) const { this->apply(dx,dy); }

    /** similar to applyDeriv */
    void applyAdjDeriv(const Vector<Scalar> & x, 
		       const Vector<Scalar> & dy,
		       Vector<Scalar> & dx) const { this->applyAdj(dy,dx); }


  public:

    LinearOp() {}
    LinearOp(const LinearOp<Scalar> & Op) {}
    virtual ~LinearOp() {}

    // getDomain and getRange inherited from Operator

    /** This function assigns to \f$ y \f$ the value \f$ A
        x\f$. Output vector y may not be aliased with input vector
        x. Applies standard sanity test, then delegates to protected
        apply method. */

    void applyOp(Vector<Scalar> const & x,
		 Vector<Scalar> & y) const {
      try {

	if (x.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: LinearOp::applyOp - input not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: LinearOp::applyOp - output not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->apply(x,y);

      }

      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyOp\n";
	throw e;
      }
    }

    /** This function assigns to \f$ y \f$ the value \f$ A^* x\f$.
	Output vector y may not be aliased with input vector
	x. Applies standard sanity test, then delegates to protected
	applyAdj method. */
   
    void applyAdjOp(Vector<Scalar> const & x,
		    Vector<Scalar> & y) const {

      try {

	if (x.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: LinearOp::applyAdjOp - input not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: LinearOp::applyAdjOp - output not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->applyAdj(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyAdjOp\n";
	throw e;
      }
    }

    /** Apply and linear combination: \f$ y = \alpha A x + \beta y \f$.
	Obvious default implementation provided for convenience, can be
	overridden to fuse loops for efficiency, if desired.  Output
	vector y may not be aliased with input vector x.
    */
    virtual void applyOp(Scalar alpha, Vector<Scalar> const & x,
			 Scalar beta, Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(this->getRange());
	// refer to public rather than deprecated protected method
	this->applyOp(x,z);
	y.linComb(alpha,z,beta);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyOp (with lincomb) \n";
	throw e;
      }
    }
    
    /** Adjoint apply and linear combination: \f$ y = \alpha A^* x + \beta y\f$. 
	Obvious default implementation provided for convenience, can
	be overridden to fuse loops for efficiency, if desired.
	Output vector y may not be aliased with input vector x.
    */
    virtual void applyAdjOp(Scalar alpha, Vector<Scalar> const & x,
			    Scalar beta, Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(this->getDomain());
	// refer to public rather than deprecated protected method
	this->applyAdjOp(x,z);
	y.linComb(alpha,z,beta);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyAdjOp (with lincomb) \n";
	throw e;
      }
    }

  };

  /** Standard construction of a LinearOp, given spaces for domain and
      range and FunctionObjects implementing the forward and adjoint
      apply... operations. */
  template<class Scalar>
  class LinearOpFO: public LinearOp<Scalar> {

  private: 

    FunctionObject & fwdfo;
    FunctionObject & adjfo;

    const Space<Scalar> & dom;
    const Space<Scalar> & rng;

    LinearOpFO();

  protected:

    LinearOp<Scalar> * clone() const {
      return new LinearOpFO<Scalar>(*this); 
    }

    void apply(const Vector<Scalar> & Input, 
	       Vector<Scalar> & Output) const {
      try {
	Output.eval(fwdfo,Input);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOpFO::applyOp\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & Input,
		  Vector<Scalar> & Output) const {
      try {
	Output.eval(adjfo,Input);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOpFO::applyAdjOp\n";
	throw e;
      }
    }

  public:

    LinearOpFO(const Space<Scalar> & _dom,
	       const Space<Scalar> & _rng,
	       FunctionObject & _fwdfo,
	       FunctionObject & _adjfo)
      : LinearOp<Scalar>(), 
	fwdfo(_fwdfo), adjfo(_adjfo), dom(_dom), rng(_rng)  {}
    LinearOpFO( const LinearOpFO<Scalar> & l)
      : LinearOp<Scalar>(l), 
	fwdfo(l.fwdfo), adjfo(l.adjfo), dom(l.dom), rng(l.rng) {}
    virtual ~LinearOpFO() {}

    virtual const Space<Scalar> & getDomain() const { return dom; }
    virtual const Space<Scalar> & getRange() const { return rng; }

    virtual ostream & write(ostream & str) const {
      str<<"Linear Operator calling function objects for image methods\n";
      str<<"domain:\n";
      dom.write(str);
      str<<"range:\n";
      rng.write(str);
      str<<"Function object implementing forward map:\n";
      fwdfo.write(str);
      str<<"Function object implementing adjoint map:\n";
      adjfo.write(str);
      return str;
    }
  };
  
  /** Invertible is a mixin interface for operators which can 
      compute inverses..
  */
  template<class Scalar>
  class Invertible {

  protected:

    /** Evaluation of linear operator inverse on constant input vector x,
	output written on mutable vector y. Accessed only through
	public applyOp method, so creator of subclasses may assume
	that input is in range, output is in domain */

    virtual void applyInv(const Vector<Scalar> & x, 
			  Vector<Scalar> & y) const = 0;

    /** Evaluation of adjoint linear operator inverse on constant input vector
	x, output written on mutable vector y. Accessed only through
	public applyOp method, so creator of subclasses may assume
	that input is in range, output is in domain */
    virtual void applyInvAdj(const Vector<Scalar> & x,
			     Vector<Scalar> & y) const = 0;

  public:
    Invertible() {}
    Invertible(const Invertible<Scalar> & Op) {}
    virtual ~Invertible() {}
  };

  /** Linear operator with inverse mapping supplied as a class method. */

  template<class Scalar>
  class LinearOpWithInverse : public LinearOp<Scalar>, public Invertible<Scalar> {
  public:
    LinearOpWithInverse() {}
    ~LinearOpWithInverse() {}

    /** This function assigns to y the value \f$ A^{-1} * x\f$, and
	is the same as solving the system \f$ A * y = x \f$ for y.
	Output vector y may not be aliased with input vector x. Applies
	standard sanity test, then delegates to protected apply
	method. */

    void applyInvOp(Vector<Scalar> const & x,
		    Vector<Scalar> & y) const {
      try {

	if (x.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: Invertible::applyInvOp - input not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: Invertible::applyInvOp - output not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->applyInv(x,y);

      }

      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyOp\n";
	throw e;
      }
    }

    /** This function assigns to y the value \f$ A^{-T} * x\f$, and
	is the same as solving the system \f$ A^{T} * y = x \f$ for y.
	Output vector y may not be aliased with input vector x. Applies
	standard sanity test, then delegates to protected applyAdj
	method. */
   
    void applyInvAdjOp(Vector<Scalar> const & x,
		       Vector<Scalar> & y) const {

      try {

	if (x.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: Invertible::applyInvAdjOp - input not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: Invertible::applyInvAdjOp - output not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->applyInvAdj(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearOp::applyAdjOp\n";
	throw e;
      }
    }

  };

  /** AdjLinearOp creates the adjoint operator of an
      LinearOp as a linear operator in its own right.  Note that
      this object points to the original linear operator.  The primary
      methods are precisely those of LinearOp; in particular, the
      image methods such as apply() and applyAdj() work by simply
      calling the appropriate methods of the original linear operator. */
  template <class Scalar>
  class AdjLinearOp: public LinearOp<Scalar> {
  
  private:
  
    const LinearOp<Scalar> & op;
  
    // default constructor--disabled
    AdjLinearOp();
    // copy constructor--private
    AdjLinearOp(const AdjLinearOp<Scalar> & a) 
      : op(a.op) {}
  
  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new AdjLinearOp<Scalar>(*this);
    }

    void apply(const Vector<Scalar> & x, 
		 Vector<Scalar> & y) const {
      try {
	op.applyAdjOp(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from AdjLinearOp::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x,
		  Vector<Scalar> & y) const {
      try {
	op.applyOp(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from AdjLinearOp::applyAdj\n";
	throw e;
      }
    }
  
  public:
  
    AdjLinearOp(const LinearOp<Scalar> & _op)
      : op(_op) {}

    ~AdjLinearOp() {}
  
    const Space<Scalar> & getDomain() const { return op.getRange(); } 
    const Space<Scalar> & getRange() const { return op.getDomain(); }

    ostream & write(ostream & str) const {
      str<<"Adjoint Operator - build on\n";
      op.write(str);
      return str;
    }
  };

  /** NormalLinearOp creates the normal operator \f$ A^*A \f$ of a LinearOp \f$ A \f$ as a
      linear operator in its own right.  Note that this object points
      to the original linear operator.  The primary methods are
      precisely those of LinearOp; in particular, the various image
      methods are implemented to call the corresponding methods of the
      original linear operator.
  */
  template<class Scalar>
  class NormalLinearOp: public LinearOp<Scalar> {
  
  private:
  
    const LinearOp<Scalar> & op;
  
    NormalLinearOp();
    NormalLinearOp(const NormalLinearOp<Scalar> & n)
      : op(n.op) {}
  
  protected:

    LinearOp<Scalar> * clone() const {
      return new NormalLinearOp<Scalar>(*this);
    }

    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(op.getRange());
	op.applyOp(x,z);
	op.applyAdjOp(z,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from NormalLinearOp::apply\n";
	throw e;
      }
    } 

    void applyAdj(const Vector<Scalar> & x,
		  Vector<Scalar> & y) const {
      try {
	this->applyOp(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from NormalLinearOp::applyAdj\n";
	throw e;
      }
    } 
    
  public:
  
    // constructs the object from a LinearOp
    NormalLinearOp(const LinearOp<Scalar> & _op)
      : op(_op) {}
  
    // destructor.
    ~NormalLinearOp() {}
  
    /** access to domain and range */
    const Space<Scalar> & getDomain() const { return op.getDomain(); }
    const Space<Scalar> & getRange() const { return op.getDomain(); }

    // report
    ostream & write(ostream & str) const {
      str<<"Normal Operator - build on\n";
      op.write(str);
      return str;
    }
  };

  /** ScaleOpFwd implementing the linear operator \f$ x\mapsto ax \f$
      where a is a scalar.  In addition to the LinearOp_d methods,
      this class has methods to access the scalar \f$ a \f$: the
      setScale() method can be used to modify that scalar whereas the
      getScale() method can be used to access it.
  */

  template<class Scalar>
  class ScaleOpFwd: public LinearOp<Scalar> {
  
  private:
  
    const Space<Scalar> & sp;
    Scalar mu;
  
    ScaleOpFwd();

  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new ScaleOpFwd<Scalar>(*this);
    }


    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	if (abs(mu-ScalarFieldTraits<Scalar>::One())<=numeric_limits<Scalar>::epsilon())
	  y.copy(x);
	else if (abs(mu)<=numeric_limits<Scalar>::epsilon())
	  y.zero();
	else
	  y.scale(mu,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScaleOpFwd::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x,
		  Vector<Scalar> & y) const {

      try {
	this->applyOp(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScaleOpFwd::applyAdj\n";
	throw e;
      }
    }  

  public:
  
    ScaleOpFwd(const Space<Scalar> & _sp, Scalar _mu)
      : sp(_sp), mu(_mu) {}

    ScaleOpFwd(const ScaleOpFwd<Scalar> & s)
      : sp(s.sp), mu(s.mu) {}

    ~ScaleOpFwd() {}
    const Space<Scalar> & getDomain() const { return sp; } 
    const Space<Scalar> & getRange() const { return sp; }
    Scalar getScale() const { return mu; }
    void setScale(Scalar _mu) { mu=_mu; }

    ostream & write(ostream & str) const {
      str<<"ScaleOpFwd - build on\n";
      sp.write(str);
      str<<"and scale factor\n";
      str<<mu<<endl;
      return str;
    }
  };

  /** ScaleOpInv implementing the linear operator \f$ x\mapsto
      \frac{1}{a}x \f$ where \f$ a \f$ is a scalar. It is built using
      the corresponding forward scale operator and offers the same
      methods as ScaleOpFwd does.
  */
  template<class Scalar>
  class ScaleOpInv: public LinearOp<Scalar> {
  
  private:
  
    const Space<Scalar> & sp;
    Scalar mu;
  
    ScaleOpInv();
    ScaleOpInv(const ScaleOpInv<Scalar> & s)
      : sp(s.sp), mu(s.mu) {}
  
  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new ScaleOpInv<Scalar>(*this);
    }

    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	if (fabs(mu-ScalarFieldTraits<Scalar>::One())<=numeric_limits<Scalar>::epsilon())
	  y.copy(x);
	else if (fabs(mu)<=numeric_limits<Scalar>::epsilon())
	  y.zero();
	else
	  y.scale(mu,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScaleOpInv::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x,
		  Vector<Scalar> & y) const {
      try {
	apply(x,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScaleOpInv::applyAdj\n";
	throw e;
      }
    }

  public:
  
    ScaleOpInv(const LinearOp<Scalar> & op,Scalar _mu)
      : sp(op.getDomain()) {
      Scalar one = ScalarFieldTraits<Scalar>::One();
      if (ProtectedDivision(one,_mu,mu)) {
	RVLException e;
	e<<"Error: ScaleOpInv constructor\n";
	e<<"reciprocal of mu = "<<_mu<<" caused zerodivide\n";
	throw e;
      }
    }
  
    ScaleOpInv(const ScaleOpFwd<Scalar> & op)
      : sp(op.getDomain()) {
      Scalar one = ScalarFieldTraits<Scalar>::One();
      if (ProtectedDivision(one,op.getScale(),mu)) {
	RVLException e;
	e<<"Error: ScaleOpInv constructor\n";
	e<<"based on ScaleOpFwd:\n";
	op.write(e);
	e<<"reciprocal of mu = "<<op.getScale()<<" caused zerodivide\n";
	throw e;
      }
    }
  
    ~ScaleOpInv() {}
  
    const Space<Scalar> & getDomain() const { return sp; } 
    const Space<Scalar> & getRange() const { return sp; }
  
    Scalar getScale() const { return mu; }
    void setScale(Scalar _mu) { mu=_mu; }

    ostream & write(ostream & str) const {
      str<<"Scale Operator - build on\n";
      sp.write(str);
      str<<"and scale factor\n";
      str<<mu<<endl;
      return str;
    }
  };

  /** LinCombLinearOp is a concrete class implementing a linear
      combination with weights \f$ w_1, w_2 \f$ of two linear
      operators \f$ Op_1 \f$ and \f$ Op_2 \f$, that is, \f$w_1 Op_1 +
      w_2 Op_2\f$

      The constructor checks that the two argument operators 
      share domain and range spaces, then creates the requested
      linear combination. Combinations of more operators can be 
      built by recursively combining pairs.

  */

  template<class Scalar>
  class LinCombLinearOp: public LinearOp<Scalar> {
  
  private:

    Scalar w1, w2;
    LinearOp<Scalar> const & op1;
    LinearOp<Scalar> const & op2;
  
    LinCombLinearOp();

  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new LinCombLinearOp<Scalar>(*this);
    }
      
    void apply(const Vector<Scalar> & x, 
		 Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(op1.getRange());
	op1.applyOp(x,z);
	op2.applyOp(x,y);
	y.linComb(w1,z,w2);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombLinearOp::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x, 
		    Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(op1.getDomain());
	op1.applyAdjOp(x,z);
	op2.applyAdjOp(x,y);
	y.linComb(w1,z,w2);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombLinearOp::applyAdj\n";
	throw e;
      }
    }
  
  public:
  
    LinCombLinearOp(Scalar _w1, LinearOp<Scalar> const & _op1,
		    Scalar _w2, LinearOp<Scalar> const & _op2)
      : w1(_w1), w2(_w2), op1(_op1), op2(_op2) {
      try {
	if ((op2.getDomain() != op1.getDomain()) ||
	    (op2.getRange() != op1.getRange())) {
	  RVLException e;
	  e<<"Error: LinCombLinearOp constructor\n";
	  e<<"incompatible domains or ranges\n";
	  e<<"\n";
	  e<<"**** operator 1:\n";
	  op1.write(e);
	  e<<"\n";
	  e<<"**** operator 2:\n";
	  op2.write(e);
	  e<<"\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombLinearOp::LinCombLinearOp\n";
	throw e;
      }
    }

    LinCombLinearOp(LinCombLinearOp const & op) 
      : w1(op.w1), w2(op.w2), op1(op.op1), op2(op.op2) {}

    ~LinCombLinearOp() {}

    const Space<Scalar> & getDomain() const {
      return op1.getDomain();
    }

    const Space<Scalar> & getRange() const {
      return op1.getRange();
    }

    ostream & write(ostream & str) const {
      str<<"LinCombLinearOp - linear combinations\n";
      str<<"\n";
      str<<"**** operator 1:\n";
      op1.write(str);
      str<<"     weight 1 = "<<w1<<"\n";
      str<<"\nfollowed by\n";
      str<<"**** operator 2:\n";
      op2.write(str);
      str<<"     weight 2 = "<<w2<<"\n";
      str<<"\n";
      return str;
    }
  };

  /** Composition of linear operators \f$ Op_1, Op_2 \mapsto Op_2
      \circ Op_1 \f$ (so subscripts indicate order of evaluation -
      that's how the constructor is organized). Alignment of domains
      and ranges checked as part of construction. */

  template<class Scalar>
  class CompLinearOp: public LinearOp<Scalar> {
  
  private:
  
    LinearOp<Scalar> const & op1;
    LinearOp<Scalar> const & op2;
  
    CompLinearOp();

  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new CompLinearOp<Scalar>(*this);
    }
      
    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(op1.getRange());
	op1.applyOp(x,z);
	op2.applyOp(z,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x, 
		    Vector<Scalar> & y) const {
      try {
	Vector<Scalar> z(op1.getRange());
	op2.applyAdjOp(x,z);
	op1.applyAdjOp(z,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::applyAdj\n";
	throw e;
      }
    }

  public:
  
    CompLinearOp(LinearOp<Scalar> const & _op1,
		 LinearOp<Scalar> const & _op2)
      : op1(_op1), op2(_op2) {
      try {
	if (op2.getDomain() != op1.getRange()) {
	  RVLException e;
	  e<<"Error: CompLinearOp constructor\n";
	  e<<"incompatible domains or ranges\n";
	  e<<"\n";
	  e<<"**** operator 1:\n";
	  op1.write(e);
	  e<<"\n";
	  e<<"**** operator 2:\n";
	  op2.write(e);
	  e<<"\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::CompLinearOp\n";
	throw e;
      }
    }

    CompLinearOp(CompLinearOp const & op) 
      : op1(op.op1), op2(op.op2) {}

    ~CompLinearOp() {}

    const Space<Scalar> & getDomain() const {
      return op1.getDomain();
    }

    const Space<Scalar> & getRange() const {
      return op2.getRange();
    }
  
    ostream & write(ostream & str) const {
      str<<"Comp LinearOperator - composition of\n";
      str<<"\n";
      str<<"**** operator 1:\n";
      op1.write(str);
      str<<"\nfollowed by\n";
      str<<"**** operator 2:\n";
      op2.write(str);
      str<<"\n";
      return str;
    }
  };

  /** for the moment, a standalone class. Symmetric merely in the sense
      that only one of two possible adjoints is supplied, the other one
      being presumed to be the same. Symmetry should be tested for any
      child class. Could eventually be made a child of OpProdDom, just
      as LinearOp is now an op.
  */
  template<typename Scalar>
  class SymmetricBilinearOp: public Writeable {

  protected:

    virtual void apply(const Vector<Scalar> & x0, 
		       const Vector<Scalar> & x1,
		       Vector<Scalar> & y) const = 0;

    virtual void applyAdj(const Vector<Scalar> & x0, 
			  const Vector<Scalar> & y,
			  Vector<Scalar> & x1) const = 0;
    
    virtual SymmetricBilinearOp<Scalar> * clone() const = 0;

  public:

#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

    virtual ~SymmetricBilinearOp() {}

    virtual Space<Scalar> const & getDomain() const = 0; 
    virtual Space<Scalar> const & getRange() const = 0; 

    /** This function assigns to \f$ y \f$ the value \f$
        A(x_0,x_1)\f$. Asserted to be same as value of \f$
        A(x_1,x_0)\f$. Output vector \f$ y \f$ may not be aliased with
        input vectors \f$ x_0, x_1 \f$. Applies standard sanity test,
        then delegates to protected apply method. */

    void applyOp(Vector<Scalar> const & x0,
		 Vector<Scalar> const & x1,
		 Vector<Scalar> & y) const {
      try {

	if (x0.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: SymmetricBilinearOp::applyOp - input 0 not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x0.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (x1.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: SymmetricBilinearOp::applyOp - input 1 not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x1.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: SymmetricBilinearOp::applyOp - output not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->apply(x0,x1,y);

      }

      catch (RVLException & e) {
	e<<"\ncalled from SymmetricBilinearOp::applyOp\n";
	throw e;
      }
    }
    
    /** This function assigns to \f$ x_1 \f$ the value \f$
        A^*(x_0,y)\f$, defined by

	\f\[ 
	\langle x_1, A^*(x_0,y) \rangle_X = \langle A(x_0,x_1), y \rangle_Y
	\f\]

	Since \f$ A \f$ is asserted to be symmetric, this definition
        of the adjoint coincides with the other obvious definition.
        Applies standard sanity test, then delegates to protected
        applyAdj method. */

    void applyAdjOp(Vector<Scalar> const & x0,
		    Vector<Scalar> const & y,
		    Vector<Scalar> & x1) const {
      try {

	if (x0.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: SymmetricBilinearOp::applyOp - input 0 not in domain\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** input space:   ***\n";
	  e<<"**********************\n";
	  x0.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (x1.getSpace() != this->getDomain()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: SymmetricBilinearOp::applyOp - output not in domain \n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** domain space:  ***\n";
	  e<<"**********************\n";
	  this->getDomain().write(e);
	  e<<"**********************\n";
	  e<<"*** output space:  ***\n";
	  e<<"**********************\n";
	  x1.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	if (y.getSpace() != this->getRange()) {
	  RVLException e;
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  e<<"Error: LinearOp::applyOp - intput 1 not in range\n";
	  e<<"**********************\n";
	  e<<"*** this operator: ***\n";
	  e<<"**********************\n";
	  this->write(e);
	  e<<"**********************\n";
	  e<<"*** range space:   ***\n";
	  e<<"**********************\n";
	  this->getRange().write(e);
	  e<<"**********************\n";
	  e<<"*** input 1 space: ***\n";
	  e<<"**********************\n";
	  y.getSpace().write(e);
	  e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	  throw e;
	}

	this->applyAdj(x0,y,x1);

      }

      catch (RVLException & e) {
	e<<"\ncalled from SymmetricBilinearOp::applyAdjOp\n";
	throw e;
      }
    }
  };

  /** LinearOp crreated by fixing the first argument in a bilinear op. */
  template<class Scalar>
  class LinearBilinearOp: public LinearOp<Scalar> {

  private:
    SymmetricBilinearOp<Scalar> const & blop;
    Vector<Scalar> const & x0;

#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

  protected:

    void apply(Vector<Scalar> const & x1,
	       Vector<Scalar> & y) const {
      try {
	blop.apply(x0,x1,y);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearBilinearOp::apply\n";
	throw e;
      }
    }

    void applyAdj(Vector<Scalar> const & y,
		  Vector<Scalar> & x1) const {
      try {
	blop.applyAdj(x0,y,x1);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinearBilinearOp::applyAdj\n";
	throw e;
      }
    }

    Operator<Scalar> * clone() const {
      return new LinearBilinearOp<Scalar>(*this);
    }

  public:

    LinearBilinearOp(SymmetricBilinearOp<Scalar> const & _blop,
		     Vector<Scalar> const & _x0) 
      : blop(_blop), x0(_x0) {}
    
    LinearBilinearOp(LinearBilinearOp<Scalar> const & lbl) 
      : blop(lbl.blop), x0(lbl.x0) {}

    Space<Scalar> const & getDomain() const { return blop.getDomain(); }
    Space<Scalar> const & getRange() const { return blop.getRange(); }

    ostream & write(ostream & str) const {
      str<<"Linear Operator wrapper around Bilinear Operator\n";
      str<<"  by fixing first argument\n";
      str<<"Bilinear Operator:\n";
      blop.write(str);
      str<<"First argument (vector):\n";
      x0.write(str);
      return str;
    }
  };
}

#endif







