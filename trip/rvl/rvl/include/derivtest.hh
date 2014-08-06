/*************************************************************************

Copyright Rice University, 2004, 2005, 2006.
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

#ifndef __RVL_DERIVTEST
#define __RVL_DERIVTEST

#include "op.hh"

namespace RVL {

  /** tests accuracy of directional derivative computation by
      comparison with centered finite difference approximation. Prints
      several columns of output, including estimated convergence rate,
      which should approach 2. Since no unversal scale rule is
      possible, users will have to adjust the limits hmin and hmax to
      suit the application.
      
      Arguments: 
      @param op operator whose deriv comp is to be tested
      @param y base point at which derivative accuracy is tested 
      @param dy direction of differentiation 
      @param str output stream on which to print
      report (use cout for terminal output) 
      @param n number of (evenly spaced) offsets h at which to
      compute divided differences. Default is 10.
      @param hmin minimum offset (default 0.1)
	@param hmax maximum offset (defaults 1.0, but these defaults
	are only useful if the operator
	is rather smooth in the intrinsic scale of \f$y\f$).
        @minrate rate of convergence lower bound - test returns true
	if the last convergence rate computed is lower than this one
	(i.e. "yes there is a problem"), else false ("no problem").

	Denoting the operator by \f$F\f$, the test computes \f$DF(y)dy\f$, and 
	compares it to the divided difference 
	\f$G(h) = (F(y+h*dy)-F(y-h*dy))/(2*h)\f$, where 
	\f$h = hmin + i*dh, dh = (hmax-hmin)/(n-1)\f$ for \f$i=0,...,n-1\f$.
	The estimated convergence rate of the divided differences to
	the derivative is 
	\f$\log(\|G(h+dh)-(h+dh)*DF(y)dy\|/\|G(h)-h*DF(y)dy\|)/
	\log((h+dh)/h)\f$, which should approach 2 as \f$h \rightarrow 0\f$.
  */

  template<typename Scalar>
  bool DerivTest(Operator<Scalar> const & op,
		 Vector<Scalar> const & y,
		 Vector<Scalar> const & p,
		 ostream & str,
		 int n = 10,
		 typename ScalarFieldTraits<Scalar>::AbsType hmin = 0.1, 
		 typename ScalarFieldTraits<Scalar>::AbsType hmax = 1.0,
		 typename ScalarFieldTraits<Scalar>::AbsType minrat = 1.95) {
    
    try {
      
      if( !y.inSpace(op.getDomain()) ) {
	RVLException e;
	e<<"Error: Operator::checkDeriv: \n";
	e<<"base vector is not in Domain\n";
	throw e;
      }
      if( !p.inSpace(op.getDomain()) ) {
	RVLException e;
	e<<"Error: Operator::CheckDeriv: \n";
	e<<"direction vector is not in Domain\n";
	throw e;
      }

      if( hmax <= hmin ) {
	typename ScalarFieldTraits<Scalar>::AbsType temp = hmax;
	hmax = hmin;
	hmin = temp;
      }

      if( hmin <= 0.0 ) {
	hmin = 0.1;
	hmax = 1.0;
      }
      if( n <= 0 ) n = 10;


      Scalar hlimit1;
      Scalar hlimit2;

      hlimit1 = op.getMaxStep( y,p );
      {
	Vector<Scalar> ptemp(op.getDomain());
	ptemp.scale( -1.0,p);
	hlimit2 = op.getMaxStep( y,ptemp );
      }

      typename ScalarFieldTraits<Scalar>::AbsType hlimit = min( abs(hlimit1),abs(hlimit2) );
      if( hlimit <= 0.0 ) {
	RVLException e;
	e<<"Error: Operator::CheckDeriv: direction is not \n";
	e<<"feasible\n";
	throw e;
      }

      if( hmax >= hlimit ) {
	hmax = 0.99*hlimit;
	if( hmin >= hmax ) hmin = hmax/n;
      }

      Vector<Scalar> x1(op.getDomain());
      Vector<Scalar> x2(op.getDomain());
      Vector<Scalar> g1(op.getRange());
      Vector<Scalar> g2(op.getRange());
      Vector<Scalar> dg(op.getRange());

      OperatorEvaluation<Scalar> Fy(op,y);
      Fy.getDeriv().applyOp(p,dg);
    
      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;

      int oldprecision = str.precision( nd );

      typename ScalarFieldTraits<Scalar>::AbsType dgnorm = dg.norm();
      int rflag = 1;

      if( dgnorm < numeric_limits<typename ScalarFieldTraits<Scalar>::AbsType >::epsilon() ) {
	rflag = 0;
	str << "DerivTest: norm of first variation = "<<dgnorm<<" is too "
	    << endl << "small; displaying absolute error" << endl;
      }
    
      str<<endl<<"Operator::checkDeriv"<<endl<<endl;
   
      if( rflag )
	str << setw(8) << "h" << setw(nd+7) << " norm of diff." << setw(nd+8) 
	    << "rel. error" << setw(nd+6)   << "convg. rate" << endl;
      else
	str << setw(8) << "h" << setw(nd+7) << " norm of diff." << setw(nd+6)
	    << "convg. rate" << endl;
      int i;
      typename ScalarFieldTraits<Scalar>::AbsType hstep = (hmax-hmin)/(n-1);
      typename ScalarFieldTraits<Scalar>::AbsType n1p=ScalarFieldTraits<Scalar>::AbsOne();
      typename ScalarFieldTraits<Scalar>::AbsType tn;
      typename ScalarFieldTraits<Scalar>::AbsType tmp =ScalarFieldTraits<Scalar>::AbsZero();

      OperatorEvaluation<Scalar> opeval1(op,x1);
      OperatorEvaluation<Scalar> opeval2(op,x2);

      for( i=n-1;i>=0;i-- )  {
	typename ScalarFieldTraits<Scalar>::AbsType h = hmin+i*hstep;
	x1.copy(y);
	x1.linComb(-h,p);
	// it is essential to use evals here - otherwise the 
	// internal workspace created for x1 won't be reinitialized
	// for evaluating on x2, as the operator clone implied by
	// the eval constructor will not take place - you can't
	// just use apply, even though it's available as a class
	// method!!!
	x2.copy(y);
	x2.linComb(h,p);
	g2.copy(opeval2.getValue());
	Scalar one = ScalarFieldTraits<Scalar>::One();
	g2.linComb(-one,opeval1.getValue());
	g2.linComb(-one,dg,one/(2.0*h));
	typename ScalarFieldTraits<Scalar>::AbsType n1 = g2.norm();
	if (i<n-1 && !ProtectedDivision<typename ScalarFieldTraits<Scalar>::AbsType>(n1p,n1,tn)) {
	  tmp = log(tn)/log((h+hstep)/h);
	  if( rflag )
	    str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
		<< setw(nd+6) << n1 << setw(nd+8) << n1/dgnorm << setw(nd+6)
		<< tmp << endl;
	  else
	    str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
		<< setw(nd+6) << n1 << setw(nd+6)
		<< tmp  << endl;
	}
	else {
	  char * jnk = new char[nd];
	  for (int ii=0;ii<nd;ii++) {
	    jnk[ii]='-';
	  }
	  jnk[nd-1]='\0';
	  if( rflag )
	    str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
		<< setw(nd+6) << n1 << setw(nd+8) << n1/dgnorm << setw(nd+6)
		<< jnk << endl;
	  else
	    str << setprecision(6) << setw(8) << h << " " << setprecision(nd)
		<< setw(nd+6) << n1 << " " << setw(nd+6)
		<< jnk  << endl;
	  delete [] jnk;
	}
	n1p=n1;
      }
      str.precision( oldprecision );
      // return value pred on last value of log ratio tmp
      if (tmp > minrat) return true;
      return false;
    }
    catch (RVLException & e) {
      e<<"\ncalled from Operator::checkDeriv\n";
      throw e;
    }

  }

}

#endif
