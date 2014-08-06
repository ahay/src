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

#ifndef __RVL_TOOLS_POLYOP_H_
#define __RVL_TOOLS_POLYOP_H_

#include "functions.hh"
#include "op.hh"

namespace RVL {

  /** This operator takes a vector of coefficients and computes the polynomial
      mapping described by these coefficients.  The coefficients are
      in ascending order of power (\f$ coef[0] + coef[1]*x + coef[2]*x^2 + \ldots\f$.
      Powers are computed elementwise, and the Vector::linComb method is used for
      addition and scaling.
  */
  template<class Scalar>
  class PolynomialOperator : public OperatorWithInvertibleDeriv<Scalar> {

  protected:
    Space<Scalar> & spc;
    std::valarray<Scalar> coef;

    /** 1-jet at a point, accessible only through OperatorEvaluation
	and subclasses.
    */
  
    /** \f$y = F(x)\f$ */
    virtual void apply(const Vector<Scalar> & x, 
		       Vector<Scalar> & y) const {
      Vector<Scalar> temp(x);
      ElementwiseMultiply<Scalar> dotstar;
      int size = coef.size();
      if( size > 0) {
	RVLAssignConst<Scalar> constantterm(coef[0]);
	y.eval(constantterm);
      }
      if( size > 1 )
	y.linComb(coef[1], temp);
      for( int i = 2; i < size; i++) {
	temp.eval(dotstar, temp, x);
	y.linComb(coef[i], temp);
      }
    }
  
    /** \f$dy = DF(x)dx\f$ */
    virtual void applyDeriv(const Vector<Scalar> & x, 
			    const Vector<Scalar> & dx,
			    Vector<Scalar> & dy) const {
      Vector<Scalar> temp(x);
      ElementwiseMultiply<Scalar> dotstar;
      int size = coef.size();
      if( size > 1) {
	RVLAssignConst<Scalar> constantterm(coef[1]);
	dy.eval(constantterm);
      } else {
	dy.zero();
      }  
      if( size > 2 )
	dy.linComb(coef[2]*Scalar(2), temp);
      for( int i = 3; i < size; i++) {
	temp.eval(dotstar, temp, x);
	dy.linComb(coef[i]*Scalar(i), temp);
      } 
      // dy has DF(x), so we need to do the multiply now.
      dy.eval(dotstar, dy,dx);
    }
  
    /** \f$dx = DF(x)^*dy\f$ */
    virtual void applyAdjDeriv(const Vector<Scalar> & x, 
			       const Vector<Scalar> & dy,
			       Vector<Scalar> & dx) const {
      applyDeriv(x,dy,dx);
    }

    /** Since this is a diagonal operator, the inverse of the derivative
	amounts to performing elementwise division.
    */
    void applyInverseDeriv(const Vector<Scalar> & x,
			   const Vector<Scalar> & dy,
			   Vector<Scalar> & dx) const {
      Vector<Scalar> temp(x);
      ElementwiseMultiply<Scalar> dotstar;
      int size = coef.size();
      if( size > 1) {
	RVLAssignConst<Scalar> constantterm(coef[1]);
	dx.eval(constantterm);
      } else {
	dx.zero();
      }  
      if( size > 2 )
	dx.linComb(coef[2]*Scalar(2), temp);
      for( int i = 3; i < size; i++) {
	temp.eval(dotstar, temp, x);
	dx.linComb(coef[i]*Scalar(i), temp);
      } 
      // dy has DF(x), so we need to do the division now.
      ElementwiseDivision<Scalar> dotslash;
      dx.eval(dotslash, dy, dx);
    }

    /** Diagonal op => symetric */
    void applyAdjInverseDeriv(const Vector<Scalar> & x,
			      const Vector<Scalar> & dx,
			      Vector<Scalar> & dy) const {
      applyInverseDeriv(x,dx,dy);
    }


    /** virtual copy contructor, also accessible only through
	OperatorEvaluation. Usually implemented with operator new and
	copy constructor of concrete child class. */
  
    virtual Operator<Scalar> * clone() const { return new PolynomialOperator(*this); }

  public:
  
    PolynomialOperator(const std::valarray<Scalar> & _coef, Space<Scalar> & _spc)
      : spc(_spc), coef(_coef) {}
    PolynomialOperator(const PolynomialOperator<Scalar> & s): spc(s.spc), coef(s.coef) {}
    ~PolynomialOperator() {}
  
    /** access to domain, range */
    virtual const Space<Scalar> & getDomain() const { return spc; }
    virtual const Space<Scalar> & getRange() const { return spc; }

    virtual void write(RVLException & e) const { 
      e << "PolynomialOperator with coefficients";
      for(int i = 0; i < coef.size(); i++) {
	e << "c[" << i << "] = " << coef[i] << '\n';
      }
    }

    virtual ostream & write(ostream & str) const { 
      str << "PolynomialOperator with coefficients";
      for(int i = 0; i < coef.size(); i++) {
	str << "c[" << i << "] = " << coef[i] << '\n';
      }
    }

  
  };


}

#endif
