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

#ifndef __RVL_FUNC
#define __RVL_FUNC

#include "op.hh"

namespace RVL {
 
  template<class Scalar>
  class FunctionalEvaluation;

  template<class Scalar>
  class Functional;

  template<class Scalar>
  class FcnlOpComp;

  template<class Scalar>
  class LinCombFunctional;

  /** Interface for scalar-valued vector functions. Provides first and
      second derivatives (gradient and Hessian) - higher derivatives
      may be added later.

      This interface has virtually no public members, other than those
      which idenify its domain. The action happens in the protected
      member functions. These get used by the Evaluation classes to
      initialize the various values (of the function itself, its
      derivatives, etc.) and any intermediate data on which these
      might depend. The clone method is also protected, and is used by
      the Evaluation constructor to create a completely independent
      copy.

      So concrete instances of this type should be written in "use-once"
      style: any and all internal data should be written once and treated
      as read-only thereafter. The only access to those functions which
      might change its internal state (all protected) is through Evaluation,
      and therefore tied to a particular evaluation point.

      Since Evaluation objects clone these, efficient implementations
      will typically allocate as much internal storage as possiblE
      dynamically, and only at the point of use. For example an array
      foo of length n used in the apply
      method would be initialized in the constructor as a pointer to
      NULL then reinitialized at the point of use:
      const void apply(...) { 
        ...  
	if (!foo) foo = new Scalar[n]; 
	...  
      } 

      A unit test for validity of the gradient computation is supplied
      as part of RVL, in the form of a standalone function
      (GradientTest). It is HIGHLY RECOMMENDED that every concrete
      Functional implementation be subjected to this test.

      Finally, this version of Function presumes that the values
      produced are of the same type as the input, that is, the
      template parameter Scalar describes the scalar field of both the
      input vector and the output value. The class is not a suitable
      abstraction for real valued functions of complex vector
      variables, for instance.

  */

  template<class Scalar>
  class Functional: public Writeable {

    friend class FunctionalEvaluation<Scalar>;
    friend class FcnlOpComp<Scalar>;
    friend class LinCombFunctional<Scalar>;

  protected:

    /** 2-jet at a point, accessible only through FunctionalEvaluation
     */
    /** \f$val = F(x)\f$ */
    virtual void apply(const Vector<Scalar> & x, 
		       Scalar & val) const = 0;

    /** \f$g = grad F(x)\f$ */
    virtual void applyGradient(const Vector<Scalar> & x, 
			       Vector<Scalar> & g) const = 0;

    /** \f$dy = Hess F(x) dx\f$ */
    virtual void applyHessian(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx, 
			      Vector<Scalar> & dy) const = 0;

    /** The export-apply methods make the protected apply methods 
	of any Functional subclass instance available to other 
	Functional subclass instances. */

    void export_apply(Functional<Scalar> const & f,
		      const Vector<Scalar> & x,
		      Scalar & val) const {
      f.apply(x,val);
    }

    void export_applyGradient(Functional<Scalar> const & f,
			      const Vector<Scalar> & x,
			      Vector<Scalar> & g) const {
      f.applyGradient(x,g);
    }

    void export_applyHessian(Functional<Scalar> const & f,
			     const Vector<Scalar> & x,
			     const Vector<Scalar> & dx, 
			     Vector<Scalar> & dy) const {
      f.applyHessian(x,dx,dy);
    }

    /** operator new - not available to general public, but
	available to children who will use it to define a 
	clone method, for example 
	 
	Version 1.0: user control
    */
#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

    /** virtual copy constructor: make a complete new copy including
	internal workspace. Usually implemented with operator new and
	copy constructor of concrete child class.
    */

    virtual Functional<Scalar> * clone() const = 0;

    void export_clone(Functional<Scalar> const & fref,
		      Functional<Scalar> ** f) const {
      try {
	if (*f) {
	  RVLException e;
	  e<<"Error: Functional::export_clone\n";
	  e<<"cannot clone to a non-null pointer\n";
	  throw e;
	}
	*f = fref.clone();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Functional::export_clone\n";
	throw e;
      }
    }
      
  public:

    Functional() {}

    Functional(const Functional<Scalar> &) {}

    virtual ~Functional() {}

    // access to domain
    virtual const Space<Scalar> & getDomain() const = 0;

    /** getMaxStep() computes the largest value of $a$ such that
	$x+a*dir$ lies in the domain of definition of the
	functional.  By default (i.e. unless this virtual
	function is implemented in a derived class), the domain is
	assumed to be the whole space, and MaxStep return the largest
	floating point number.  This method is used in some of the
	optimization algorithms and provides a partial solution to the
	problem that many functionals are not defined on an entire
	vector space, but rather only on a subset.
	It can be overridden in derived classes */
    virtual Scalar getMaxStep(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx) const {
      return numeric_limits<Scalar>::max();
    }

    /** Formerly built-in gradient accuracy test. Now delegates to
	standalone function.

    bool checkGradient(const Vector<Scalar> & y,
		       const Vector<Scalar> & p,
		       ostream & str,
		       int n=10,
		       Scalar hmin=0.1, 
		       Scalar hmax=1.0,
		       Scalar minrate = 0.0) {
      return GradientTest(*this,y,p,str,n,hmin,hmax,minrate);
    }
    */
    /** Formerly built-in hessian accuracy test. Now delegates to
	standalone function
    bool checkHessian(const Vector<Scalar> & y,
		      const Vector<Scalar> & p,
		      ostream & str,
		      int n=10,
		      Scalar hmin=0.1,
		      Scalar hmax=1.0) {
      return HessianTest(*this,y,p,str,n,hmin,hmax);
    }

    . */
  };

  template<class Scalar>
  class FunctionalProductDomainEvaluation;

  /** A specialization of Functional which also has partial derivatives. */
  template<class Scalar>
  class FunctionalProductDomain: public Functional<Scalar> {

    friend class FunctionalProductDomainEvaluation<Scalar>;

  protected:

    /** applyPartialGradient() computes the component of the gradient
	corresponding to the ith component of the independent variable.
	This is a virtual function and must be implemented in any derived
	class. */
    virtual void applyPartialGradient(int i,
				      const Vector<Scalar> & x,
				      Vector<Scalar> & g) const = 0;

    /** applyGradient() computes the gradient of the function. This method is
	implemented in this (the base) class, using PartialGradient(),
	and can be (but need not be) re-implemented in a derived class. */
    virtual void applyGradient(const Vector<Scalar> & x,
			       Vector<Scalar> & g) const {
      try {
	Components<Scalar> cx(x);
	applyPartialGradient(0,cx[0],g);
	Vector<Scalar> tmp(g);
	for (int i=1;i<cx.getSize();i++) {
	  applyPartialGradient(i,cx[i],tmp);
	  g.linComb(1.0,tmp);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalProductDomain::applyGradient\n";
	throw e;
      }
    }

    /** applyHessianBlock() computes the image of the (i,j) block of the
	Hessian on derivative on dxi, giving dxj.  If this method is not
	implemented, then the default HessianBlock() method can be
	used; the block of the Hessian will then be an instance of 
	FunctionalDefaultHessianBlock. */
    virtual void applyHessianBlock(int i,
				   int j,
				   const Vector<Scalar> & x,
				   const Vector<Scalar> & dxi,
				   Vector<Scalar> & dxj) const = 0;

    /** applyHessian() computes the image of the Hessian of the
	function. This method is implemented in this (the base) class,
	using applyHessianBlock, but may be re-implemented in a
	derived class. NB: as of 15.09.06, NOT IMPLEMENTED.*/ 
    virtual void applyHessian(const Vector<Scalar> & x,
			      const Vector<Scalar> & yin,
			      Vector<Scalar> & yout) const {
      RVLException e;
      e<<"Error: FunctionalProductDomain::applyHessian\n";
      e<<"this method not yet defined - complain to management!\n";
      throw e;
    }
    virtual FunctionalProductDomain<Scalar> * clonePD() const = 0;
    Functional<Scalar> * clone() const { return clonePD(); }

  public:

    // Default constructor.
    FunctionalProductDomain() {}
    // Copy constructor.
    FunctionalProductDomain(const FunctionalProductDomain<Scalar> &) {}
    // Destructor.
    virtual ~FunctionalProductDomain() {}

    /** getDomain() returns a reference to the domain of the underlying
	function. */
    const Space<Scalar> & getDomain() const { 
      return getProductDomain(); 
    }
    /** ProductDomain() returns a reference to the domain of the
	underlying function, as an instance of HCL_ProductSpace.
	This method is implemented using Domain(), and need not be
	re-implemented in a derived class. */
    virtual const ProductSpace<Scalar> & getProductDomain() const = 0;

    // This method checks the analytic second derivatives against the
    // second order finite-difference formula
    //
    //     \[D_{ij}f(x)(\delta x_j,\delta x_i) \doteq 
    //          f(x+\tilde{\delta x_i}+\tilde{\delta x_j})+
    //          f(x-\tilde{\delta x_i}-\tilde{\delta x_j})-
    //          f(x+\tilde{\delta x_i})-
    //          f(x-\tilde{\delta x_i})-
    //          f(x+\tilde{\delta x_j})-
    //          f(x-\tilde{\delta x_j})+2f(x)\].
    //
    //  Here $\tilde{\delta x_i}$ is the vector in the domain of $f$ with every
    //  component equal to zero except the $i$th, which equals $\delta x_i$.
  
    /** Formerly built-in test for accuracy of block Hessian
	computation - general case. Now deferred to standalone
	function.
    bool checkHessianBlock(const Vector<Scalar> & y,
			   const Vector<Scalar> & pi,
			   const Vector<Scalar> & pj,
			   int i,
			   int j,
			   ostream & str,
			   int n=10,
			   Scalar hmin=0.1,
			   Scalar hmax=1.0) {
      return HessianBlockTestNonDiag(*this,y,pi,pj,i,j,str,n,hmin,hmax);
    }
 */
    /** Formerly built-in test for accuracy of block Hessian
	computation - diagonal case. Now deferred to standalone
	function.
    bool checkHessianBlock(const Vector<Scalar> & y,
			   const Vector<Scalar> & pi,
			   int i,
			   ostream & str,
			   int n=10,
			   Scalar hmin=0.1,
			   Scalar hmax=1.0) {
      return HessianBlockTestDiag(*this,y,pi,i,n,hmin,hmax);
    }
 */
  };

  template<class Scalar>
  class HessianEvaluation;

  template<class Scalar>
  class HessianBlockEvaluation;

  /** Evaluation is a pair of a (clone of a) Functional and an
      evaluation point Vector, stored by reference. The Functional is
      re-cloned automatically when the evaluation point Vector changes
      internal state.
  */
  template<class Scalar>
  class FunctionalEvaluation: public Writeable {

    friend class HessianEvaluation<Scalar>;
    friend class HessianBlockEvaluation<Scalar>;
    friend class FcnlOpComp<Scalar>;
    typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  private:

    // Vector and functional are independent copies
    const Functional<Scalar> & fref;
    WatchedVecRef<Scalar> wx;
    mutable Functional<Scalar> * f;

    mutable Scalar val;
    mutable bool applied;
    mutable Vector<Scalar> grad;
    mutable bool gapplied;
    mutable NormRetType gnorm;
    mutable bool gnormapplied;

    Components<Scalar> cg;
    HessianEvaluation<Scalar> hess;

    // disabled
    FunctionalEvaluation();

    void reset() const {
      try {
	if (f) delete f;
	f = fref.clone();
	applied=false;
	gapplied=false;
	gnormapplied=false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::reset\n";
	throw e;
      }
    }

  protected:

    // access through HessianEvaluation
    void applyHessian(const Vector<Scalar> & yin,
		      Vector<Scalar> & yout) const {
      try {
	if (wx.update()) reset();
	f->applyHessian(wx.get(),yin,yout);
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::applyHessian()\n";
	throw e;
      }
    }
    
    /** The next two functions throw exceptions if the 
	referenced functional is not a FcnlProdDom. They are 
	accessed only by HessianBlockEval, via FcnlProdDomEval,
	which provides compile time type-safety in addition
	to the run-time type checking built into these methods. 
    */
    
    const ProductSpace<Scalar> & getProductDomain() {
      
      try {
	const FunctionalProductDomain<Scalar> & pf =
	  dynamic_cast<const FunctionalProductDomain<Scalar> &>(fref);
	return pf.getProductDomain();
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: FunctionalEvaluation::getProductDomain\n";
	e<<"referenced Functional does not have ProductSpace domain\n";
	throw e;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getProductDomain\n";
	throw e;
      }
    }

    const Vector<Scalar> & getGradientBlock(int i) {
      try {
	if (wx.update()) reset();
	if (!gapplied) {
	  grad.relLock();
	  f->applyGradient(wx.get(),grad);
	  gapplied=true;
	  grad.setLock();
	}
	return cg[i];
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getGradientBlock\n";
	throw e;
      }
    }

    /** applyHessianBlock() computes the image of the (i,j) block of the
	Hessian on derivative on dxi, giving dxj.
    */
    void applyHessianBlock(int i,
			   int j,
			   const Vector<Scalar> & dxi,
			   Vector<Scalar> & dxj) const {
      try {
	if (wx.update()) reset();
	FunctionalProductDomain<Scalar> * pf = NULL;
	if ((pf = dynamic_cast<FunctionalProductDomain<Scalar> *>(f))) {
	  pf->applyHessianBlock(i,j,wx.get(),dxi,dxj);
	}
	else {
	  RVLException e;
	  e<<"Error: FunctionalEvaluation::applyHessianBlock\n";
	  e<<"referenced Functional does not have ProductSpace domain\n";
	  e<<"so Hessian block structure not defined\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation";
	e<<"::applyHessianBloack()\n";
	throw e;
      }
    }

    /** operator new - not available to general public, but
	available to children and friends 

	Version 1.0: user control
    */

#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

  public:

    /** main constructor - combines const references to Functional,
	Vector. The Functional is cloned, and the clone is used as
	working copy. The const reference is copied (along with that
	to the evaluation point Vector) and used whenever fresh
	Functional clones are needed. */
    FunctionalEvaluation(const Functional<Scalar> & _f, 
			 const Vector<Scalar> & x)
      : fref(_f), wx(x), f(_f.clone()), 
	applied(false), grad(fref.getDomain()), 
	gapplied(false), gnormapplied(false),
	cg(grad), hess(*this)
	 {
      grad.zero();
      
      if (x.getSpace() != fref.getDomain()) {
	RVLException e;
	e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	e<<"Error: FunctionalEvaluation constructor\n";
	e<<"-- input not in domain of Functional \n";
	e<<"**********************\n";
	e<<"*** this fnctnl:   ***\n";
	e<<"**********************\n";
	fref.write(e);
	e<<"**********************\n";
	e<<"*** domain space:  ***\n";
	e<<"**********************\n";
	fref.getDomain().write(e);
	e<<"**********************\n";
	e<<"*** input space:   ***\n";
	e<<"**********************\n";
	x.getSpace().write(e);
	e<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	throw e;
      }

    }

    FunctionalEvaluation(const FunctionalEvaluation<Scalar> & ev)
      : wx(ev.wx), f(ev.fref.clone()), fref(ev.fref), 
	applied(false),
	grad(fref.getDomain()),
	cg(grad), gapplied(false), gnormapplied(false),
	hess(*this) { grad.zero(); }

    virtual ~FunctionalEvaluation() { 
      try { 
	if (f) delete f;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionEvaluation destructor\n";
	throw e;
      }
    }

    /** access to domain */
    const Space<Scalar> & getDomain() const {
      try { return fref.getDomain(); }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getDomain\n";
	throw e;
      }
    }

    /** Returns the maximum scalar alpha for which
	x+alpha*dx is still in the domain of the functional.
    */
    Scalar getMaxStep(const Vector<Scalar> & dx) const {
      const Vector<Scalar> & fruit = getPoint();
      return fref.getMaxStep(fruit,dx);
    }

    /** const reference to evaluation point */
    Vector<Scalar> & getPoint() { return wx.get(); }
    Vector<Scalar> const & getPoint() const { return wx.get(); }

    /** extract value of functional at evaluation point. Checks to see
	if latter has been updated; if so, updates internal copy of
	value. */
    Scalar getValue() const {
      try {
	if (wx.update()) {
	  reset();
	}
	if (!applied) {
	  f->apply(wx.get(),val);
	  applied = true;
	}
	return val;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getValue()\n";
	throw e;
      }
    }

    /** returns a reference to the gradient of the functional at the
	current point. WARNING: It is not safe to save a reference to
	the gradient. This method should be called to ensure proper
	recalculation of the gradient when the point changes.
    */
    Vector<Scalar> const & getGradient() const {
      try {
	if (wx.update()) {
	  reset();
	}
	if (!gapplied) {
	  grad.zero();
	  f->applyGradient(wx.get(),grad);
	  gapplied=true;
	}
	return grad;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getGradient\n";
	throw e;
      }
    }
    
    /** Returns the norm of the gradient of the functional at the current
	point.  WARNING:  It is not safe to save a reference to the gradient.
	This method (or the getGradient() method) should be called to ensure 
	proper recalculation of the gradient when the point changes.
    */
    Scalar getGradientNorm() const {
      try {
	if (wx.update()) {
	  reset();
	}
	if (!gnormapplied) {
	  gnorm = getGradient().norm();
	  gnormapplied = true;
	}
	return gnorm;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getGradientNorm\n";
	throw e;
      }
    }

    /** Returns a LinearOp representing the Hessian of the functional.
	This LinearOp is a lightweight object that refers back to the
	evaluation, and get destroyed in the evaluation destructor.
	It is safe to store a reference to the Hessian for the lifetime of
	the evaluation.
    */
    LinearOp<Scalar> const & getHessian() const {
      try { return hess; }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalEvaluation::getHessian\n";
	throw e;
      }
    }
    
    /** provided to enable extraction of subclass attributes via a cast, 
	from the current internal copy of the underlying Functional. 
    */
    Functional<Scalar> const & getFunctional() const { return *f; }

    ostream & write(ostream & str) const {
      str<<"Functional Evaluation:"<<endl;
      str<<"  functional"<<endl;
      fref.write(str);
      str<<"  evaluated at\n";
      wx.get().write(str);
      return str;
    }
  };

  /** The Hessian Evaluation is a lightweight implementation of the 
      LinearOp interface which refers back to a FunctionalEvaluation
      to implement all methods.
  */
  template<class Scalar>
  class HessianEvaluation: public LinearOp<Scalar> {

  private:
    
    FunctionalEvaluation<Scalar> & fx;
    
    // disabled
    HessianEvaluation();
    HessianEvaluation(const HessianEvaluation<Scalar> & h)
      : fx(h.fx) {}

  protected:
    
    LinearOp<Scalar> * clone() const {
      return new HessianEvaluation<Scalar>(*this);
    }

    // image, application, MatVec product, whatever
    void apply(const Vector<Scalar> & yin, 
		 Vector<Scalar> & yout) const {
      try {
	fx.applyHessian(yin,yout);
      }
      catch (RVLException & e) {
	e<<"\ncalled from HessianEvaluation::apply()\n";
	throw e;
      }
    }
    
    // image of adjoint (transpose) operator
    void applyAdj(const Vector<Scalar> & yin,
		    Vector<Scalar> & yout) const {
      try {
	fx.applyHessian(yin,yout);
      }
      catch (RVLException & e) {
	e<<"\ncalled from HessianEvaluation::applyAdj()\n";
	throw e;
      }
    }
    
  public:

    HessianEvaluation(FunctionalEvaluation<Scalar> & _fx): fx(_fx) {}
    ~HessianEvaluation() {}

    // access to domain and range
    const Space<Scalar> & getDomain() const { return fx.getDomain(); }
    const Space<Scalar> & getRange() const { return fx.getDomain(); }

    /** report to stream */
    ostream & write(ostream & str) const {
      str<<"Hessian operator"<<endl;
      str<<"part of functional evaluation"<<endl;
      fx.write(str);
      return str;
    }
  };

  /** A specialization of FunctionalEvaluation which accesses the
      additional partial derivatives of the FunctionalProductDomain
      class. NOTE 15.09.06: this class is UNDER CONSTRUCTION and
      should be used WITH CAUTION. */
  template<class Scalar>
  class FunctionalProductDomainEvaluation: 
    public FunctionalEvaluation<Scalar> {

  private:

    HessianBlockEvaluation<Scalar> deriv;

    // disabled
    FunctionalProductDomainEvaluation();
    FunctionalProductDomainEvaluation
    (const FunctionalProductDomainEvaluation<Scalar> &);
  
  public:
  
    FunctionalProductDomainEvaluation(FunctionalProductDomain<Scalar> & _f, 
				      const Vector<Scalar> & _x)
      : FunctionalEvaluation<Scalar>(_f,_x), deriv(*this) {}

    ~FunctionalProductDomainEvaluation() {}

    const Vector<Scalar> & getGradientBlock(int i) {
      try {
	return FunctionalEvaluation<Scalar>::getGradientBlock(i);
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalProductDomainEvaluation";
	e<<"::getGradientBlock\n";
	throw e;
      }
    }

    /** This seems like a bad implementation.  At the moment, only
	one block may be accessed at a time, and references to previous
	blocks in fact all point to the same object.  Use with caution.
    */
    const LinearOp<Scalar> & getHessianBlock(int i, int j) { 
      try {
	deriv.setBlock(i,j);
	return deriv;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FunctionalProductDomainEvaluation";
	e<<"::getHessianBlock\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"Functional Evaluation with Product Domain; as"<<"\n";
      return FunctionalEvaluation<Scalar>::write(str);
    }
  };

  /** Accesses a single block of the Hessian of a functional over a
      product domain.  A single object represents every block of the
      entire Hessian, as the block access indices can be reset by
      using the setBlock() method.  This object has a lifetime
      strictly controlled by the functional with which it is
      associated and cannot survive independently.

      ADP: Is there someway to restrict the use of this to be only
      available inside the FPDEval?
  */
  template<class Scalar>
  class HessianBlockEvaluation: public LinearOp<Scalar> {
    
  private:

    FunctionalProductDomainEvaluation<Scalar> & fx;
    int i;
    int j;

    // disabled
    HessianBlockEvaluation();
    HessianBlockEvaluation(const HessianBlockEvaluation<Scalar> & h)
      : fx(h.fx), i(h.i), j(h.j) {}

  protected:

    HessianBlockEvaluation(FunctionalProductDomainEvaluation<Scalar> & _fx)
      : fx(_fx), i(0), j(0) {}
    
    LinearOp<Scalar> * clone() const {
      return new HessianBlockEvaluation<Scalar>(*this);
    }

    // image, application, MatVec product, whatever
    void apply(const Vector<Scalar> & y, 
	       Vector<Scalar> & z) const {
      try {
	fx.applyHessianBlock(i,j,y,z);
      }
      catch (RVLException & e) {
	e<<"\ncalled in HessianBlockEvaluation::apply\n";
	throw e;
      }
    }

    // image of adjoint (transpose) operator
    void applyAdj(const Vector<Scalar> & y,
		    Vector<Scalar> & z) const {
      try {
	fx.applyHessianBlock(i,j,y,z);
      }
      catch (RVLException & e) {
	e<<"\ncalled in HessianBlockEvaluation::applyAdj\n";
	throw e;
      }
    }

  public:

    ~HessianBlockEvaluation() {}

    const Space<Scalar> & getDomain() const { 
      try {
	return (fx.getProductDomain())[i];
      }
      catch (RVLException & e) {
	e<<"\ncalled from HessianBlockEvaluation::getDomain()\n";
	throw e;
      }
    }

    const Space<Scalar> & getRange() const { 
      try {
	return (fx.getProductDomain())[j];
      }
      catch (RVLException & e) {
	e<<"\ncalled from HessianBlockEvaluation::getDomain()\n";
	throw e;
      }
    }

    void setBlock(int ii, int jj) { i=ii; j=jj; }

    ostream & write(ostream & str) const {
      str<<"Hessian block operator"<<"\n";
      str<<"block ("<<i<<", "<<j<<")\n";
      str<<"part of functional evaluation"<<"\n";
      fx.write(str);
      return str;
    }
  };

  /** LinCombFunctional is a concrete class implementing a linear
      combination of two or more Functional instances.

      To construct the Functional M = 2L1-3.7L2+9.0L3, where
      L1, L2 and L3 are Functional instances,

      LinCombFunctional<float> M;
      M.setNext(2.0,L1);
      M.setNext(-3.7,L2);
      M.setNext(9.0,L3);

      Because linear combinations of two Functionals occur frequently,
      and for backwards compatibility, a special constructor is
      provided for this case:

      LinCombFunctional<float> M(2.0,L1,-3.7,L2);

      constructs M = 2L1-3.7L2.

      Post-construction initialization permits linear combinations of
      arbitrary length, but requires that the user signal when the
      initialization is complete, so that well-defined behaviour
      results. As in OpComp, this class uses the convention that the
      first non-initialization method call flags the completion of the
      linear combination. No subsequent modifications are permitted.

      The domain is determined by the Functional in the first
      summand. All subsequent summands are checked for matching domain
      as they are submitted.
  */

  template<class Scalar>
  class LinCombFunctional: public Functional<Scalar> {
  private:
  
    mutable std::vector<Functional<Scalar> *> fnvec;
    mutable std::vector<Scalar> wtvec;
    mutable bool applied;

    /** Run-time initialization */
    void setNext(Scalar a, Functional<Scalar> & fn) {
      try {
	if (applied) {
	  RVLException e;
	  e<<"Error: LinCombFunctional::setNext\n";
	  e<<"object already initialized - non-initialization method called\n";
	  e<<"further alteration to object data not allowed\n";
	  throw e;
	}
	if (fnvec.size() > 0) {
	  if (fn.getDomain() != fnvec[0]->getDomain()) {
	    RVLException e;
	    e<<"Error: LinCombOp::setNext\n";
	    e<<"domain of input Functional incompatible with reference (summand 0)\n";
	    e<<"*** input Functional:\n";
	    fn.write(e);
	    e<<"*** reference Functional:\n";
	    fnvec[0]->write(e);
	    throw e;
	  }
	}
	wtvec.push_back(a);
	fnvec.push_back(fn.clone());
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional::setNext\n";
	throw e;
      }
    }

  protected:

    void apply(const Vector<Scalar> & x, 
	       Scalar & val) const {
      try {
   
	if (fnvec.size()<1) {
	  RVLException e;
	  e<<"Error: LinCombFcnl::apply\n";
	  e<<"not initialized\n";
	  throw e;
	}
	applied = true;
	export_apply(*(fnvec[0]),x,val);
	val *= wtvec[0];
	if (fnvec.size() > 1) {
	  Scalar tmp;
	  for (int i=1; i<fnvec.size(); i++) {
	    export_apply(*(fnvec[i]),x,tmp);
	    val += wtvec[i]*tmp;
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional::apply\n";
	throw e;
      }
    }

    void applyGradient(const Vector<Scalar> & x,
		       Vector<Scalar> & g) const {
      try {
	if (fnvec.size()<1) {
	  RVLException e;
	  e<<"Error: LinCombFcnl::applyGradient\n";
	  e<<"not initialized\n";
	  throw e;
	}
	applied = true;
	export_applyGradient(*(fnvec[0]),x,g);
	g.scale(wtvec[0]);
	if (fnvec.size() > 1) {
	  Vector<Scalar> tmp(fnvec[0]->getDomain());
	  for (int i=1; i<fnvec.size(); i++) {
	    export_applyGradient(*(fnvec[i]),x,tmp);
	    g.linComb(wtvec[i],tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional::applyDeriv\n";
	throw e;
      }
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx,
		      Vector<Scalar> & dy) const {
      try {
	if (fnvec.size()<1) {
	  RVLException e;
	  e<<"Error: LinCombOp::applyHessian\n";
	  e<<"not initialized\n";
	  throw e;
	}
	applied = true;
	export_applyHessian(*(fnvec[0]),x,dx,dy);
	dy.scale(wtvec[0]);
	if (fnvec.size() > 1) {
	  Vector<Scalar> tmp(fnvec[0]->getDomain());
	  for (int i=1; i<fnvec.size(); i++) {
	    export_applyHessian(*(fnvec[i]),x,dx,tmp);
	    dy.linComb(wtvec[i],tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional::applyHessian\n";
	throw e;
      }
    }

    Functional<Scalar> * clone() const { 
      applied = true;
      return new LinCombFunctional<Scalar>(*this); 
    }

  public:
  
    LinCombFunctional() {}
    LinCombFunctional(LinCombFunctional<Scalar> const & fn) {
      try {
	for (int i=0;i<fn.fnvec.size(); i++) {
	  fnvec.push_back(fn.fnvec[i]->clone());
	  wtvec.push_back(fn.wtvec[i]);
	}	
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional copy constructor\n";
	throw e;
      }
    }

    LinCombFunctional(Scalar a1,const Functional<Scalar> & fn1,
		      Scalar a2,const Functional<Scalar> & fn2)
      : fnvec(2), wtvec(2) {
      try {
	if (fn1.getDomain() != fn2.getDomain()) {      
	  RVLException e;
	  e<<"Error: LinCombFunctional pair constructor\n";
	  e<<"domains do not matcth\n";
	  e<<"first Functional:\n";
	  fn1.write(e);
	  e<<"second Functional:\n";
	  fn2.write(e);
	  throw e;
	}
	fnvec[0] = fn1.clone();
	fnvec[1] = fn2.clone();
	wtvec[0] = a1;
	wtvec[1] = a2;
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional pair constructor\n";
	throw e;
      }
    }

    ~LinCombFunctional() {
      for (int i=0;i<fnvec.size(); i++) if (fnvec[i]) delete fnvec[i];
    }

    /** access to domain and range */
    const Space<Scalar> & getDomain() const {
      try {
	if (fnvec.size()<1) {
	  RVLException e;
	  e<<"Error: LinCombFunctional::getDomain\n";
	  e<<"object not initialized\n";
	  throw e;
	}
	applied = true;
	return fnvec[0]->getDomain();
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinCombFunctional::getDomain\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"LinCombFunctional: linear combination of scalar-valued functions\n";
      if (fnvec.size()<1) {
	str<<"not initialized\n";
      }
      else {
	for (int i=0;i<fnvec.size();i++) {
	  str<<" --- Functional "<<i<<" with weight "<<wtvec[i]<<"\n";
	  fnvec[i]->write(str);
	}
      }
      return str;
    }
  };

  /**  This class implements the Functional interface by combining the
       operation of three FunctionObjects and a domain Space.

       Note that FunctionObject instances have persistent state and
       thus are not required to yield the same result when applied to
       the same data. Thus only FunctionObjects which do not change
       results on repeated evaluation at the same instance are usable
       in this construction.

       Storage of intermediate results likely introduces
       inefficiencies. The use of this class is primarily for rapid
       prototyping or the implementation of functionals with no
       intermediate data between functional, gradient, and hessian
       evaluations. It is possible to get around this limitation by
       having these three FOs reference a common external object.

       The functional evaluation is performed by a
       FunctionObjectScalarRedn which returns a Scalar value.  The
       gradient computation is performed by a binary FunctionObject,
       which will be given the gradient (target) vector as the first
       input and the point x (data) vector as the second input.  The
       hessian computation is performed by a ternary FunctionObject
       which will be given the target vector dy first, followed by the
       point x and direction dx \f$dy = H(x)*dx\f$.
  */

  template<class Scalar, class DataType = Scalar>
  class StdFOFunctional : public Functional<Scalar> {
  protected:
    FunctionObjectScalarRedn<Scalar> & f;
    FunctionObject & gradf;
    FunctionObject & hessf;
    Space<Scalar> const & dom;    

    /** 2-jet at a point, accessible only through FunctionalEvaluation
     */
    /** \f$val = F(x)\f$ */
    virtual void apply(const Vector<Scalar> & x, 
		       Scalar & val) const {
      //replace with component-wise statements (D.S. 04.10.11)
          x.eval(f);
          val = f.getValue(); 
	  //Components<Scalar> cx(x);
	  //int csize = cx.getSize();
	  //val = ScalarFieldTraits<Scalar>::Zero();
	  //for (int i=0; i< csize; i++) {
	  //cx[i].eval(f);
	  //val += f.getValue();
	  //}
    }

    /** \f$g = grad F(x)\f$ */
    virtual void applyGradient(const Vector<Scalar> & x, 
			       Vector<Scalar> & g) const {
      g.eval(gradf,x);
    }

    /** \f$dy = Hess F(x) dx\f$ */
    virtual void applyHessian(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx, 
			      Vector<Scalar> & dy) const {
      dy.eval(hessf,x,dx);
    }

    /** virtual copy constructor: make a complete new copy including
	internal workspace. Virtual to permit
	override in child class.
    */

    virtual Functional<Scalar> * clone() const {
      return new StdFOFunctional<Scalar,DataType>(*this);
    }

    StdFOFunctional();
  public:
    
    /** main constructor - takes references to functionals defining
	main protected methods, domain space. */
    StdFOFunctional( FunctionObjectScalarRedn<Scalar> & f_,
		     FunctionObject & gradf_,
		     FunctionObject & hessf_,
		     const Space<Scalar> & dom_)
      : f(f_), gradf(gradf_), hessf(hessf_), dom(dom_) {}

    StdFOFunctional( const StdFOFunctional<Scalar, DataType> & s)
      : f(s.f), gradf(s.gradf), hessf(s.hessf), dom(s.dom) {}
    
    ~StdFOFunctional() {}

    // access to domain
    virtual const Space<Scalar> & getDomain() const { return dom; }

    virtual ostream & write(ostream & str) const {
      str << "StdFOFunctional with f = ";
      f.write(str);
      str << "\n gradf = ";
      gradf.write(str);
      str << "\n hessf = ";
      hessf.write(str);
      return str;
    }
  };

  /*
  template<class Scalar>
  bool Functional<Scalar>::checkHessian(const Vector<Scalar> & y,
					const Vector<Scalar> & p,
					ostream & str,
					int n,
					Scalar hmin,
					Scalar hmax) {
    try {
      if (!y.inSpace(getDomain())) {
	RVLException e; e<<"Error in Functional::checkGradient: \n";
	e<<"base vector is not in Domain\n";
	throw e;
      }
      if (!p.inSpace(getDomain())) {
	RVLException e; e<<"Error in Functional::checkGradient: \n";
	e<<"direction vector is not in Domain\n";
	throw e;
      }
   
      if (hmax <= hmin) {
	Scalar temp = hmax;
	hmax = hmin;
	hmin = temp;
      }
      if (hmin <= 0.0) {
	hmin = 0.1;
	hmax = 1.0;
      }
      if (n <= 0)
	n = 10;
      
      Scalar hlimit1;
      Scalar hlimit2;
      hlimit1 = getMaxStep(y,p);
      {
	Vector<Scalar> ptemp(getDomain());
	ptemp.scale(-1.0,p);
	hlimit2 = getMaxStep(y,ptemp);
      }
      Scalar hlimit = min(hlimit1,hlimit2);
      if (hlimit <= 0.0) {
	RVLException e; e<<"Error in Functional::checkHessian: \n";
	e<<"direction is not feasible\n";
	throw e;
      }
      if (hmax >= hlimit) {
	hmax = 0.99*hlimit;
	if (hmin >= hmax)
	  hmin = hmax/n;
      }

      Vector<Scalar> x1(getDomain());
      Vector<Scalar> x2(getDomain());
      Vector<Scalar> g1(getDomain());
      Vector<Scalar> g2(getDomain());
      Vector<Scalar> dg(getDomain());

      FunctionalEvaluation<Scalar> Fy(*this,y);
      Fy.getHessian().applyOp(p,dg);

      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;

      int oldprecision = str.precision(nd);

      Scalar dgnorm = dg.norm();
      int rflag = 1;
      if (dgnorm < numeric_limits<Scalar>::epsilon()) {
	rflag = 0;
	str << "Functional::checkHessian: norm of gradient is too "
	    << endl << "small; displaying absolute error" << endl;
      }
      
      str<<endl<<"Hessian Computation Check"<<endl<<endl;
      
      if (rflag)
	str << setw(8) << "h" << setw(nd+7) << "norm of diff."
	    << setw(nd+8) << "norm" << setw(nd+6) << "Rel. Err."
	    << endl;
      else
	str << setw(8) << "h" << setw(nd+7) << "norm of diff."
	    << setw(nd+6) << "norm" << endl;
      int i;
      Scalar hstep = (hmax-hmin)/(n-1);
      for (i=n-1;i>=0;i--) {
	Scalar h = hmin+i*hstep;
	x1.copy(y);
	x1.linComb(-h,p);
	FunctionalEvaluation<Scalar> F1(*this,x1);
	x2.copy(y);
	x2.linComb(h,p);
	FunctionalEvaluation<Scalar> F2(*this,x2);
	g2.copy(F2.getGradient());
	g2.linComb(-1.0,F1.getGradient());
	g2.scale(1.0/(2.0*h));
	g2.linComb(-1.0,dg);
	Scalar n1 = g2.norm();
	if (rflag)
	  str << setprecision(6) << setw(8) << h << " "
	      << setprecision(nd) << setw(nd+6) << n1 << " " 
	      << setw(nd+6) << dgnorm << " " << setw(nd+6)
	      << n1/dgnorm << endl;
	else
	  str << setprecision(6) << setw(8) << h << " "
	      << setprecision(nd) << setw(nd+6) << n1 << " " 
	      << setw(nd+6) << dgnorm << " " << endl;
      }
      str.precision(oldprecision);
      return 0;
    }
    catch (RVLException & e) {
      e<<"\ncalled from Functional::checkHessian\n";
      throw e;
    }
  }

  template<class Scalar>
  bool FunctionalProductDomain<Scalar>::checkHessianBlock
  (const Vector<Scalar> & y,
   const Vector<Scalar> & pi,
   const Vector<Scalar> & pj,
   int i,
   int j,
   ostream & str,
   int n,
   Scalar hmin,
   Scalar hmax) {
    try {
      if (i==j)
	return checkHessianBlock(y,pi,i,str,n,hmin,hmax);
      
      if (i < 1 || i > getProductDomain().getSize() || j < 1 ||
	  j > getProductDomain().getSize()) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"invalid index\n";
	throw e;
      }
      if (!y.inSpace(getDomain())) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"base vector is not in Domain\n";
	throw e;
      }
      if (!pi.inSpace(getProductDomain()[i])) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"direction vector pi is not in Domain\n";
	throw e;
      }
      if (!pj.inSpace(getProductDomain()[j])) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"direction vector pj is not in Domain\n";
	throw e;
      }

      if (hmax <= hmin) {
	Scalar temp = hmax;
	hmax = hmin;
	hmin = temp;
      }
      if (hmin <= 0.0) {
	hmin = 0.1;
	hmax = 1.0;
      }
      if (n <= 0)
	n = 10;
      
      // There are now six directions to check: pi+pj,-pi-pj,pi,-pi,pj,-pj
      Scalar hlimit;
      {
	Vector<Scalar> yp(getDomain());
	yp.zero();
	Components<Scalar> cy(yp);
	yp[i].copy(pi);
	yp[j].copy(pj);
	hlimit = getMaxStep(y,yp);
	yp.negate();
	hlimit = min(getMaxStep(y,yp),hlimit);
	yp.zero();
	yp[i].copy(pi);
	hlimit = min(getMaxStep(y,yp),hlimit);
	yp.negate();
	hlimit = min(getMaxStep(y,yp),hlimit);
	yp.zero();
	yp[j].copy(pj);
	hlimit = min(getMaxStep(y,yp),hlimit);
	yp.negate();
	hlimit = min(getMaxStep(y,yp),hlimit);
      }
      if (hlimit <= 0.0) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: direction is not feasible\n";
	throw e;
      }
      if (hmax >= hlimit) {
	hmax = 0.99*hlimit;
	if (hmin >= hmax)
	  hmin = hmax/n;
      }

      Vector<Scalar> x1(getDomain());
      Scalar v0,dv;

      FunctionalProductDomainEvaluation<Scalar> Fy(*this,y);
      Vector<Scalar> pi1(getProductDomain()[i]);
      Fy.getHessianBlock(i,j).applyOp(pj,pi1);
      dv = pi.inner(pi1);

      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;
      int oldprecision = str.precision(nd);

      Scalar dvmag = abs(dv);
      int rflag = 1;
      if (dvmag < numeric_limits<Scalar>::epsilon()) {
	rflag = 0;
	str << "FunctionalProductDomain::CheckHessianBlock: norm of "
	  "second variation is too "
	    << endl << "small; displaying absolute error" << endl;
      }
      if (rflag)
	str << setw(5) << "h" << setw(nd+6) << "2nd variation" << setw(nd+6)
	    << "cent. diff." << setw(nd+4) << "Rel. Err." << endl;
      else
	str << setw(5) << "h" << setw(nd+6) << "2nd variation" << setw(nd+6)
	    << "cent. diff." << setw(nd+4) << "Error" << endl;
      int ii;
      Scalar hstep = (hmax-hmin)/(n-1);
      Scalar val;
      Components<Scalar> cx1(x1);
      for(ii=n-1;ii>=0;ii--) {
	Scalar v = 2.0*Fy.getValue();
	Scalar h = hmin+ii*hstep
;
	x1.copy(y);
	cx1[i].linComb(h,pi);
	cx1[j].linComb(h,pj);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v += Fp.getValue();
	}

	x1.copy(y);
	cx1[i].linComb(-h,pi);
	cx1[j].linComb(-h,pj);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v += Fp.getValue();
	}
	
	x1.copy(y);
	cx1[i].linComb(h,pi);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v -= Fp.getValue();
	}

	x1.copy(y);
	cx1[i].linComb(-h,pi);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v -= Fp.getValue();
	}

	x1.copy(y);
	cx1[j].linComb(h,pj);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v -= Fp.getValue();
	}

	x1.copy(y);
	cx1[j].linComb(-h,pj);
	{
	  FunctionalEvaluation<Scalar> Fp(*this,x1);
	  v -= Fp.getValue();
	}

	v /= (2.0*h*h);
	Scalar n1 = abs(v-dv);
	
	if (rflag)
	  str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
	      << setw(nd+6) << dv << " " << setw(nd+6)
	      << v << " " << setw(nd+6) << n1/dvmag << endl;
	else
	  str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
	      << setw(nd+6) << dv << " " << setw(nd+6)
	      << v << " " << setw(nd+6) << n1 << endl;
      }
      str.precision(oldprecision);
      return 0;
    }
    catch (RVLException & e) {
      e<<"\ncalled from FunctionalProductDomain::checkHessianBlock\n";
      throw e;
    }
  }

  template<class Scalar>
  bool FunctionalProductDomain<Scalar>::checkHessianBlock
  (const Vector<Scalar> & y,
   const Vector<Scalar> & pi,
   int i,
   ostream & str,
   int n,
   Scalar hmin,
   Scalar hmax) {
    try {
      if (i < 1 || i > i > getProductDomain()) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"invalid index\n";
	throw e;
      }
      if (!y.inSpace(getDomain())) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"base vector is not in Domain\n";
	throw e;
      }
      if (!pi.inSpace(getProductDomain()[i])) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: \n";
	e<<"direction vector pi is not in Domain\n";
	throw e;
      }

      if (hmax <= hmin) {
	Scalar temp = hmax;
	hmax = hmin;
	hmin = temp;
      }
      if (hmin <= 0.0) {
	hmin = 0.1;
	hmax = 1.0;
      }
      if (n <= 0)
	n = 10;

      Scalar hlimit;
      {
	Vector<Scalar> yp(getDomain());
	yp.zero();
	Components<Scalar> cy(yp);
	cy[i].copy(pi);
	hlimit = getMaxStep(y,yp);
	yp.negate();
	hlimit = min(getMaxStep(y,yp),hlimit);
      }
      if (hlimit <= 0.0) {
	RVLException e; e<<"Error in FunctionalProductDomain::";
	e<<"checkHessianBlock: direction is not feasible\n";
	throw e;
      }
      if (hmax >= hlimit) {
	hmax = 0.99*hlimit;
	if (hmin >= hmax)
	  hmin = hmax/n;
      }

      Vector<Scalar> x1(getDomain());
      Scalar v0,dv;

      FunctionalEvaluation<Scalar> Fy(*this,y);
      Vector<Scalar> pi1(getProductDomain()[i]);
      Fy.getHessianBlock(i,i).applyOp(pi,pi1);
      dv = pi.inner(pi1);

      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;
      int oldprecision = str.precision(nd);

      Scalar dvmag = abs(dv);
      int rflag = 1;
      if (dvmag < numeric_limits<Scalar>::epsilon()) {
	rflag = 0;
	str << "FunctionalProductDomain_d::CheckHessianBlock: norm of "
	  "second variation is too "
	    << endl << "small; displaying absolute error" << endl;
      }

      if (rflag)
	str << setw(5) << "h" << setw(nd+6) << "2nd variation" << setw(nd+6)
	    << "cent. diff." << setw(nd+4) << "Rel. Err." << endl;
      else
	str << setw(5) << "h" << setw(nd+6) << "2nd variation" << setw(nd+6)
	    << "cent. diff." << setw(nd+4) << "Error" << endl;
      Scalar hstep = (hmax-hmin)/(n-1);
      Scalar val;
      Components<Scalar> cx1(x1);
      int ii;
      for (ii=n-1;ii>=0;ii--) {
	Scalar v = -2.0*v0;
	Scalar h = hmin+ii*hstep;
	x1.copy(y);
	cx1[i].linComb(h,pi);
	{
	  FunctionalEvaluation<Scalar> F1(*this,x1);
	  v += F1.getValue();
	}

	x1.copy(y);
	cx1[i].linComb(-h,pi);
	{
	  FunctionalEvaluation<Scalar> F1(*this,x1);
	  v += F1.getValue();
	}

	v /= (h*h);
	Scalar n1 = abs(v-dv);

	if (rflag)
	  str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
	      << setw(nd+6) << dv << " " << setw(nd+6)
	      << v << " " << setw(nd+6) << n1/dvmag << endl;
	else
	  str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
	      << setw(nd+6) << dv << " " << setw(nd+6)
	      << v << " " << setw(nd+6) << n1 << endl;
      }
      str.precision(oldprecision);
      return 0;
    }
    catch (RVLException & e) {
      e<<"\ncalled from FunctionalProductDomain::checkHessianBlock\n";
      throw e;
    }
  }

  */


  /** Null-functional.  \f$f(x) = 0 \f$.  */
  template<class Scalar>
  class NullFunctional : public Functional<Scalar> {

  private:    

    Space<Scalar> const & dom; 

    NullFunctional();

  protected:
    /** \f$val = F(x)\f$ */
    virtual void apply(const Vector<Scalar> & x, 
		       Scalar & val) const {
      val = 0;
    }
    
    /** \f$g = grad F(x)\f$ */
    virtual void applyGradient(const Vector<Scalar> & x, 
			       Vector<Scalar> & g) const {
      g.zero();
    }
    
    /** \f$dy = Hess F(x) dx\f$ */
    virtual void applyHessian(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx, 
			      Vector<Scalar> & dy) const {
      dy.zero();
    }

  public:

    NullFunctional(const Space<Scalar> & sp)
      : dom(sp) {}

    NullFunctional(const Functional<Scalar> & _f)
      : dom(_f.getDomain()) {}

    ~NullFunctional() { }

    /** virtual copy constructor: make a complete new copy including
	internal workspace. Virtual to permit override.
    */

    virtual Functional<Scalar> * clone() const {
      return new NullFunctional<Scalar>(*this);
    }

    // access to domain
    virtual const Space<Scalar> & getDomain() const {
      return dom;
    }

    virtual ostream & write(ostream & str) const {
      str << "NullFunctional on space\n ";
      dom.write(str);
      str << " \n";
      return str;
    }
  };

  /** This handle class creates the composite of a functional and an
      operator \f$f(G(x))\f$, using the protected services of Operator
      and Functional.  While it can access the protected Functional
      methods through inheritance, it must be a friend of the base
      Operator class.

      The basic methods remain virtual so that further attributes can
      be added in subclasses.

      Note that the value is available once the apply method is called
      for the first time, with additional computations. However the
      intermediate vector holding the output of the operator factor is
      not, to avoid potential heavyweight storage.
  */
  template<class Scalar>
  class FcnlOpComp: public Functional<Scalar> {

  private:

    Functional<Scalar> const & f;
    Operator<Scalar> const & op;
    mutable FunctionalEvaluation<Scalar> * fneval;
    mutable OperatorEvaluation<Scalar> * opeval;

  protected:

    void apply(const Vector<Scalar> & x, 
	       Scalar & val) const {
      try {
	if (!opeval) 
	  opeval = new OperatorEvaluation<Scalar>(op,x);
	if (!fneval)
	  fneval = new FunctionalEvaluation<Scalar>(f,opeval->getValue());
	//	cerr<<"FcnlOpComp::apply -> getValue()\n";
	val=fneval->getValue();
	// delete fneval; fneval=NULL;
	// delete opeval; opeval=NULL;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FcnlOpComp::apply\n";
	throw e;
      }
    }

    /** \f$g = grad F(x)\f$ */
    virtual void applyGradient(const Vector<Scalar> & x, 
			       Vector<Scalar> & g) const {
      try {
	if (!opeval)
	  opeval = new OperatorEvaluation<Scalar>(op,x);
	if (!fneval)
	  fneval = new FunctionalEvaluation<Scalar>(f,opeval->getValue());
	Vector<Scalar> const & gtmp = fneval->getGradient();
	opeval->getDeriv().applyAdjOp(gtmp,g);
	//	delete fneval; fneval=NULL;
	//	delete opeval; opeval=NULL;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FcnlOpComp::applyGradient\n";
	throw e;
      }
    }

    virtual void applyHessian(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx, 
			      Vector<Scalar> & dy) const {

      try {
	if (!opeval) 
	  opeval = new OperatorEvaluation<Scalar>(op,x);
	if (!fneval) 
	  fneval = new FunctionalEvaluation<Scalar>(f,opeval->getValue());
	Vector<Scalar> tmp1(op.getRange());
	Vector<Scalar> tmp2(op.getRange());
	opeval->getDeriv().applyOp(dx,tmp1);
	fneval->getHessian().applyOp(tmp1,tmp2);
        opeval->getDeriv().applyAdjOp(tmp2,dy);
	//	delete fneval; fneval=NULL;
	//	delete opeval; opeval=NULL;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FcnlOpComp::applyHessian\n";
	throw e;
      }
    }

    /** virtual copy constructor: make a complete new copy including
	internal workspace. Usually implemented with operator new and
	copy constructor of concrete child class.
    */

    virtual Functional<Scalar> * clone() const { 
      //      cerr<<"FcnlOpComp::clone\n";
      return new FcnlOpComp<Scalar>(*this);
    }

  public:

    FcnlOpComp(const Functional<Scalar> & fref,
	       const Operator<Scalar> & opref)
      : f(fref), op(opref),
	fneval(NULL), opeval(NULL) {}

    FcnlOpComp(const FcnlOpComp<Scalar> & c) 
      : f(c.f), op(c.op), fneval(NULL), opeval(NULL) {
#if 0
	  cerr<<"\n\n***************************\n";
	  cerr<<"NEW FCNLOPCOMP:\n";
	  cerr<<"FCNL:\n";
	  f.write(cerr);
	  cerr<<"\nOP:\n";
	  op.write(cerr);
	  cerr<<"\n***************************\n\n";
#endif
    }

    ~FcnlOpComp() {
      if (fneval) delete fneval;
      if (opeval) delete opeval;
    }

    /** access to domain */
    const Space<Scalar> & getDomain() const {
      try { 
	return op.getDomain(); 
      }
      catch (RVLException & e) {
	e<<"\ncalled from FcnlOpComp::getDomain\n";
	throw e;
      }
    }

    virtual Scalar getMaxStep(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx) const {
      try {
	/* this is a tricky bit - in fact it is not clear
	   how in general you should interpret this. A natural
	   method might be the next block of code, i.e. linearize
	   the correspondence in direction between dom(op) and dom(f).
	   however it's clearly not right, so for the moment the 
	   domain of f plays no role in this computation, which is 
	   not right either.
	   
	   Since it's so tricky, we've left it virtual - override 
	   with anything sensible.
	*/
	/*
	Scalar os = op.get().getMaxStep(x,dx);
	if (!opapplied) {
	  op.get().apply(x,opx);
	  opapplied=true;
	}
	Vector<Scalar> tmp1(op.get().getRange());
	op.get().applyDeriv(x,dx,tmp1);
	Scalar fs = f.get().getMaxStep(opx,tmp1);
	Scalar st = min(fs,os);
	return st;
	*/
	return op.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from FcnlOpComp::getMaxStep\n";
	throw e;
      }
    }

    /** provides access to any available attributes of latest operator 
	evaluation */
    OperatorEvaluation<Scalar> const & getOpEval() const { 
      if (opeval) return *opeval;
      RVLException e;
      e<<"Error: FcnlOpComp::getOpEval\n";
      e<<"  apparently not yet evaluated so operator evaluation component\n";
      e<<"  not available\n";
      throw e;
    }

    /** provides access to any available attributes of latest operator 
	evaluation */
    FunctionalEvaluation<Scalar> const & getFcnlEval() const { 
      if (fneval) return *fneval;
      RVLException e;
      e<<"Error: FcnlOpComp::getFcnlEval\n";
      e<<"  apparently not yet evaluated so functional evaluation component\n";
      e<<"  not available\n";
      throw e;
    }

    ostream & write(ostream & str) const {
      str<<"FcnlOpComp: Functional-Operator composition\n";
      str<<"operator:\n";
      op.write(str);
      str<<"followed by functional:\n";
      f.write(str);
      return str;
    }
    
  };

}

#endif





