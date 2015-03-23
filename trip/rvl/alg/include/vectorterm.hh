
/*************************************************************************

Copyright Rice University, 2004.
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

#ifndef __RVLALG_VECTOR_TERMINATOR
#define __RVLALG_VECTOR_TERMINATOR

#include "functional.hh"
#include "ls.hh"
#include "alg.hh"

/**  vectorterm.H
     Terminators which involve manipulating Vectors.  Includes
     using the built-in Vector operations as well as passing vectors
     to other objects.
*/

namespace RVLAlg {
  using namespace RVL;
/** Terminator which takes a unary function object,
    a vector, and a tolerance.

    STOP: \f$f(x) < tol\f$
*/

template <class Scalar>
class UnaryThresholdTerminator: public Terminator {
  
public:
  UnaryThresholdTerminator( FunctionObjectScalarRedn<Scalar>& tf, 
			    Vector<Scalar> & tx,
			    Scalar ttol) 
    : f(tf), x(tx), tol(ttol) {}

  virtual bool query() { 
    x.eval(f);
    Scalar temp = f.getValue();
    return (temp < tol);
  }

protected:
  FunctionObjectScalarRedn<Scalar> & f;
  Vector<Scalar> & x;
  Scalar tol;
  
};

/** Terminator which takes a binary function object, two
    vectors, and a tolerance.

    STOP: \f$f(x,y) < tol\f$
*/

template <class Scalar>
class BinaryThresholdTerminator: public Terminator {
  
public:
  BinaryThresholdTerminator( FunctionObjectScalarRedn<Scalar>& tf, 
			     Vector<Scalar> & tx,
			     Vector<Scalar> & ty,
			     Scalar ttol) 
    : f(tf), x(tx), y(ty), tol(ttol) {}

  virtual bool query() { 
    Scalar temp;
    x.eval(f, y);
    temp = f.getValue();
    return (temp < tol);
  }

protected:
  FunctionObjectScalarRedn<Scalar> & f;
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  Scalar tol;  
};

/** Terminator which takes a ternary function object, three
    vectors, and a tolerance.

    STOP: \f$f(x,y,z) < tol\f$
*/

template <class Scalar>
class TernaryThresholdTerminator: public Terminator {
  
public:
  TernaryThresholdTerminator( FunctionObjectScalarRedn<Scalar>& tf, 
			      Vector<Scalar> & tx,
			      Vector<Scalar> & ty,
			      Vector<Scalar> & tz,
			      Scalar ttol) 
    : f(tf), x(tx), y(ty), z(tz), tol(ttol) {}

  virtual bool query() { 
    Scalar temp;
    x.eval(f,y,z);
    temp = f.getValue();
    return (temp < tol);
  }

protected:
  FunctionObjectScalarRedn<Scalar> & f;
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  Vector<Scalar> & z;
  Scalar tol;
};

/** Terminator which takes a vector and a tolerance.

    STOP: \f$\|\|x\|\| < tol\f$
*/

template< class Scalar >
class NormThresholdTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  NormThresholdTerminator( const Vector<Scalar> & tx, NormRetType ttol): x(tx), tol(ttol) {}

  virtual bool query() {
    return (x.norm() < tol);
  }

protected:
  const Vector<Scalar> & x;
  NormRetType tol;
};

/** Terminator which takes a vector and a tolerance.

    STOP: \f$\|\|x\|\|^{2} < tol\f$

    note: this is slightly cheaper than finding ||x|| < tol 
*/

template< class Scalar >
class Norm2ThresholdTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  Norm2ThresholdTerminator( Vector<Scalar> & tx, NormRetType ttol): x(tx), tol(ttol) {}

  virtual bool query() {
    return (x.norm2() < tol);
  }

protected:
  Vector<Scalar> & x;
  NormRetType tol;
};

/** Terminator which takes two vectors and a tolerance.

    STOP: \f$\|\|x - y\|\| < tol\f$
*/


template< class Scalar >
class DiffThresholdTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  DiffThresholdTerminator( Vector<Scalar> & tx, Vector<Scalar> & ty, NormRetType ttol)
    : x(tx), y(ty), tol(ttol) {
    if( ! x.inSameSpace(y) ) {
      RVLException e; e << "Error in DiffThresholdTerminator constructor: Vectors not in same space.";
      throw e;
    }
  }

  virtual bool query() {
    Vector<Scalar> temp(x.getSpace());
    temp.lincomb(1.0, x, -1.0, y);
    return (temp.norm() < tol);
  }

protected:
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  NormRetType tol;
};


/** Terminator which takes two vectors and a tolerance.

    STOP: \f$\|\|x - y\|\|^{2} < tol\f$

    note: this is slightly cheaper than finding \f$\|\|x - y\|\| < tol\f$ 
*/


template< class Scalar >
class Diff2ThresholdTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  Diff2ThresholdTerminator( Vector<Scalar> & tx, Vector<Scalar> & ty, NormRetType ttol)
    : x(tx), y(ty), tol(ttol) {
    if( ! x.inSameSpace(y) ) {
      RVLException e; e << "Error in DiffThresholdTerminator constructor: Vectors not in same space.";
      throw e;
    }
  }

  virtual bool query() {
    Vector<Scalar> temp(x.getSpace());
    temp.lincomb(1.0, x, -1.0, y);
    return (temp.norm2() < tol);
  }

protected:
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  NormRetType tol;
};

/** Terminator which takes two vectors and a tolerance.

    STOP: \f$\left< x, y \right> < tol\f$
*/

template< class Scalar >
class IPThresholdTerminator: public Terminator {
public:
  IPThresholdTerminator( Vector<Scalar> & tx, Vector<Scalar> & ty, Scalar ttol)
    : x(tx), y(ty), tol(ttol) {
    if( ! x.inSameSpace(y) ) {
	RVLException e; e << "Error in DiffThresholdTerminator constructor: Vectors not in same space.";
      throw e;
    }
  }

  virtual bool query() {
    Scalar temp = x.inner(y);
    return (temp < tol);
  }

protected:
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  Scalar tol;
};

/** Terminator which takes two vectors and a tolerance.

    STOP: \f$\| \left< x , y \right> \| < tol\f$
 */

template< class Scalar >
class AbsIPThresholdTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  AbsIPThresholdTerminator( Vector<Scalar> & tx, Vector<Scalar> & ty, NormRetType ttol)
    : x(tx), y(ty), tol(ttol) {
    if( ! x.inSameSpace(y) ) {
      RVLException e; e << "Error in DiffThresholdTerminator constructor: Vectors not in same space.";
      throw e;
    }
  }
  ~AbsIPThresholdTerminator(){}

  virtual bool query() {
    Scalar temp = x.inner(y);
    return (fabs(temp) < tol);
  }

protected:
  Vector<Scalar> & x;
  Vector<Scalar> & y;
  NormRetType tol;
};

/** A terminator which checks for a stationary point in 
    a functional.  Stops when
\f[ \|\|\nabla f(x)\|\| <= tol \f]

    Note: For efficiency, we store \f$tol^2\f$ and compare it to
    the square of the norm, since this is ussually cheaper.
*/
template <class Scalar> 
class NormGradientTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;  

  /** Takes the current x, the functional f, and the tolerance.  
      Notice that all inputs except for the tolerance are references,
      and similar internal data members are also references.  Thus, external
      changes to these objects will affect the terminator.  This
      is the intended behavior. 
      
  */
  NormGradientTerminator( Vector<Scalar> & x, 
			  Functional<Scalar> & f, 
			  NormRetType tol) 
    : tol_(tol*tol), fx_(f,x)
  {}

  /** Returns true when norm of the gradient is less than tol */
  bool query() {
    return( fx_.getGradient().norm2() <= tol_ );
  }

protected:
  NormRetType tol_; 
  FunctionalEvaluation<Scalar>  fx_;
};


/** Terminator which takes vectors x0 and x and a max for the norm difference r. Returns false if

    \f$\|x - x0\| < r\f$

    Returns true otherwise, <i>and</i> projects x onto the ball of
    radius r centered at x0, that is, replaces x by

    \f$ x \leftarrow x0 + \frac{r}{\|x-x_0\|}(x-x_0) \f$

    In order to avoid possible false positive return due to roundoff
    when query() is called repeatedly on the same data, scale
    \f$x-x_0\f$ by slightly smaller factor than returned by floating
    point division.
 */

template< class Scalar >
class DiffBallProjTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  DiffBallProjTerminator( Vector<Scalar> const & tx, Vector<Scalar> & ty, NormRetType _maxstep, ostream & _str = cout)
    : x(tx), y(ty), maxstep(_maxstep), str(_str), res(false) {
    if( ! x.inSameSpace(y) ) {
      RVLException e; e << "Error in DiffBallProjTerminator constructor: Vectors not in same space.";
      throw e;
    }
  }

  virtual bool query() {
    Vector<Scalar> temp(x.getSpace());
    temp.copy(y);
    temp.linComb(-1.0,x);
    NormRetType tn = temp.norm();
    res = (tn > maxstep);
    if (res) {
      NormRetType scfac=ScalarFieldTraits<NormRetType>::One();
      if (ProtectedDivision(maxstep,tn,scfac)) {
	RVLException e; 
	e<<"Error: zerodivide in DiffBalProjTerminator::query\n";
	e<<"maxstep = "<<maxstep<<"\n";
	e<<"temp.norm() = "<<tn<<"\n";
	throw e;
      }
	  
      temp.scale((ScalarFieldTraits<NormRetType>::One() - 100*numeric_limits<NormRetType>::epsilon())*scfac);
      y.copy(x);
      y.linComb(1.0,temp);
    }
    return res;
  }

  /** for post facto use - access to result without recomputation */
  bool static_query() { return res; }

protected:
  Vector<Scalar> const & x;
  Vector<Scalar> & y;
  NormRetType maxstep;
  bool res;
  ostream & str;
};

/** Terminator which takes x and a max for the norm r. Returns false if

    \f$\|x\| < r\f$

    Returns true otherwise, <i>and</i> projects x onto the ball of
    radius r centered at 0, that is, replaces x by

    \f$ x \leftarrow \frac{r}{\|x\|}x \f$

    Special case of DiffBallProjTerminator with x0=0.

    In order to avoid possible false positive return due to roundoff
    when query() is called repeatedly on the same data, scale
    \f$x\f$ by slightly smaller factor than returned by floating
    point division.
 */

template< class Scalar >
class BallProjTerminator: public Terminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  BallProjTerminator(Vector<Scalar> & ty, NormRetType _maxstep,
		     ostream & _str = cout)
    : y(ty), maxstep(_maxstep), queryres(false), str(_str) {}

  virtual bool query() {
    if (queryres) return true;
    NormRetType tn = y.norm();
    bool res = (tn > maxstep);
    if (res) {
      NormRetType scfac=ScalarFieldTraits<NormRetType>::One();
      if (ProtectedDivision(maxstep,tn,scfac)) {
	RVLException e; 
	e<<"Error: zerodivide in DiffBalProjTerminator::query\n";
	e<<"maxstep = "<<maxstep<<"\n";
	e<<"temp.norm() = "<<tn<<"\n";
	throw e;
      }
      y.scale((ScalarFieldTraits<NormRetType>::One() - 100*numeric_limits<NormRetType>::epsilon())*scfac);
      str<<"RVLAlg::BallProjTerminator::query: trust region truncation applied\n";
      str<<"  untruncated solution norm = "<<tn<<" trust radius = "<<maxstep<<"\n";      
    }
    queryres=res;
    return res;
  }

protected:
  Vector<Scalar> & y;
  NormRetType maxstep;
  mutable bool queryres;
  ostream & str;
};


}

#endif
