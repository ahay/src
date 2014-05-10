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

#ifndef __RVL_GRADTEST
#define __RVL_GRADTEST

#include "functional.hh"

namespace RVL {
 
  /** This implementation of the gradient test should be sufficient
      for all Functionals. Tests n finite difference steps from y in
      direction p, ranging from hmin to hmax. Compares computed rate
      of change <grad f(y), p> with centered finite difference
      approximations (f(y+h*p)-f(y-h*p))/(2*h) for n values of h
      ranging from hmin to hmax. Computes nominal convergence rate by
      Richardson extrapolation, prints result to ostream specified in
      arg list. Convergence rate should converge to 2.0 for
      well-chosen n, range of steps h.
  */
  template<class Scalar>
  bool GradientTest(Functional<Scalar> const & f,
		    const Vector<Scalar> & y,
		    const Vector<Scalar> & p,
		    ostream & str,
		    int n=11,
		    Scalar hmin=0.1, 
		    Scalar hmax=1.0,
		    Scalar minrat=1.95) {

    try {

      bool success = false;

      if (!y.inSpace(f.getDomain())) {
	RVLException e; e<<"Error in GradientTest: \n";
	e<<"base vector is not in Domain\n";
	throw e;
      }
      if (!p.inSpace(f.getDomain())) {
	RVLException e; e<<"Error in GradientTest: \n";
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
      if (n < 2) { 
	RVLException e;
	e<<"RVL::GradientTest:\n";
	e<<"  test not performed because number of samples = "<<n<<" too small\n";
	throw e;
      }

      Scalar hlimit1;
      Scalar hlimit2;
      hlimit1 = f.getMaxStep(y,p);

      {
	Vector<Scalar> ptemp(f.getDomain());
	Scalar one = ScalarFieldTraits<Scalar>::One();
	ptemp.scale(-one,p);
	hlimit2 = f.getMaxStep(y,ptemp);
      }

      Scalar hlimit = min(hlimit1,hlimit2);
      if (hlimit <= 0.0) {
	RVLException e; e<<"Error in GradientTest:\n";
	e<<" direction is not feasible\n";
	throw e;
      }
      if (hmax >= hlimit) {
	hmax = 0.99*hlimit;
	if (hmin >= hmax) hmin = hmax/n;
      }

      Vector<Scalar> x1(f.getDomain());
      Vector<Scalar> x2(f.getDomain());
      Vector<Scalar> dx(f.getDomain());

      FunctionalEvaluation<Scalar> F1(f,x1);
      FunctionalEvaluation<Scalar> F2(f,x2);

      Scalar s = 0.0;
      {
	FunctionalEvaluation<Scalar> Fy(f,y);
        s = (Fy.getGradient()).inner(p);
      }

      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;

      int oldprecision = str.precision(nd);

      int rflag = 1;

      if (abs(s) < numeric_limits<Scalar>::epsilon()) {
	rflag = 0;
	str << "GradientTest: norm of first "
	  "variation is too small; displaying absolute error"
	    << endl;
      }

      str<<endl<<"Gradient Computation Check"<<endl<<endl;
   
      if (rflag)
	str << setw(8) << "h" << setw(nd+7) << " norm of diff."
	    << setw(nd+8) << "rel. error" << setw(nd+6)
	    << "convg. rate" << endl;
      else
	str << setw(8) << "h" << setw(nd+7) << " norm of diff."
	    << setw(nd+6) << "convg. rate" << endl;
      int i;
      Scalar hstep = (hmax-hmin)/(n-1);
      Scalar n1p=1.0;
      Scalar tn;

      for(i=n-1;i>=0;i--)  {

	Scalar h = hmin+i*hstep;
     	x1.copy(y);
	x2.copy(y);
	x1.linComb(-h,p);
	x2.linComb(h,p);
	/*
	cerr<<"\ny:\n";
	y.write(cerr);
	cerr<<"\np:\n";
	p.write(cerr);
	cerr<<"\nx1:\n";
	x1.write(cerr);
	cerr<<"\nx2:\n";
	x2.write(cerr);
	cerr<<endl;
	*/
	Scalar ds = (F2.getValue()-F1.getValue())/(2.0*h);
	// cerr<<"h = "<<h<<" val1 = "<<F1.getValue()<<" val2 = "<<F2.getValue()<<" ds = "<<ds<<endl;
	Scalar n1 = abs(s-ds);
	if (i<n-1 && !ProtectedDivision(n1p,n1,tn)) {
	  Scalar rat = log(tn)/log((h+hstep)/h);
	  if (rat > minrat) success = true;
	  if (rflag)
	    str << setprecision(6) << setw(8) << h << " "
		<< setprecision(nd) << setw(nd+6) << n1
		<< setw(nd+8) << n1/abs(s) << setw(nd+6)
		<< rat << endl;
	  else
	    str << setprecision(6) << setw(8) << h << " "
		<< setprecision(nd) << setw(nd+6) << n1
		<< setw(nd+6) << rat
		<< endl;
	}
	else {
	  char * jnk = new char[nd];
	  for (int ii=0;ii<nd;ii++) {
	    jnk[ii]='-';
	  }
	  jnk[nd-1]='\0';
	  if (rflag)
	    str << setprecision(6) << setw(8) << h << " "
		<< setprecision(nd) << setw(nd+6) << n1
		<< setw(nd+8) << n1/abs(s) << setw(nd+6)
		<< jnk << endl;
	  else
	    str << setprecision(6) << setw(8) << h << " "
		<< setprecision(nd) << setw(nd+6) << n1 << " "
		<< setw(nd+6) << jnk  << endl;
	  delete [] jnk;
	}
	n1p=n1;

      }
      str.precision(oldprecision);
      return success;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GradientTest\n";
      throw e;
    }
  }

}

#endif
