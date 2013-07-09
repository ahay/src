/*************************************************************************

Copyright Rice University, 2004-2013
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

#ifndef __RVL_SCANTEST
#define __RVL_SCANTEST

#include "functional.hh"

namespace RVL {
 
  /** Computes values of an RVL::Functional along a line segment in its 
      domain, and writes them nicely formatted to an output stream.
  */
  template<class Scalar>
  void Scan(Functional<Scalar> const & f,
	    const Vector<Scalar> & y,
	    const Vector<Scalar> & p,
	    int n=11,
	    Scalar hmin=-ScalarFieldTraits<Scalar>::One(), 
	    Scalar hmax=ScalarFieldTraits<Scalar>::One(),
	    ostream & str=cout) {

    try {

      if (!y.inSpace(f.getDomain())) {
	RVLException e; e<<"Error in Scan: \n";
	e<<"  base vector is not in Domain\n";
	throw e;
      }
      if (!p.inSpace(f.getDomain())) {
	RVLException e; e<<"Error in Scan: \n";
	e<<"  direction vector is not in Domain\n";
	throw e;
      }
      if (n<0) {
	RVLException e;
	e<<"Error in Scan\n";
	e<<"  requested negative number of steps n="<<n<<"\n";
	throw e;
      }

      if (hmax > f.getMaxStep(y,p)) {
	RVLException e;
	e<<"Error in Scan\n";
	e<<"  max step "<<hmax<<" leaves domain of input function\n";
	throw e;
      }	

      if (hmax > f.getMaxStep(y,p)) {
	RVLException e;
	e<<"Error in Scan\n";
	e<<"  max step "<<hmax<<" leaves domain of input function\n";
	throw e;
      }
      if (hmin > f.getMaxStep(y,p)) {
	RVLException e;
	e<<"Error in Scan\n";
	e<<"  min step "<<hmax<<" leaves domain of input function\n";
	throw e;
      }	
      if (hmax<0) {
	Vector<Scalar> ptemp(f.getDomain());
	Scalar one = ScalarFieldTraits<Scalar>::One();
	ptemp.scale(-one,p);
	if (-hmax > f.getMaxStep(y,ptemp)) {
	  RVLException e;
	  e<<"Error in Scan\n";
	  e<<"  max step "<<hmax<<" leaves domain of input function\n";
	  throw e;
	}
      }
      if (hmin<0) {
	Vector<Scalar> ptemp(f.getDomain());
	Scalar one = ScalarFieldTraits<Scalar>::One();
	ptemp.scale(-one,p);
	if (-hmin > f.getMaxStep(y,ptemp)) {
	  RVLException e;
	  e<<"Error in Scan\n";
	  e<<"  max step "<<hmax<<" leaves domain of input function\n";
	  throw e;
	}
      }
      Vector<Scalar> x(f.getDomain());
      FunctionalEvaluation<Scalar> feval(f,x);

      int nd;
      if (numeric_precision<Scalar>()==1) nd = 8;
      if (numeric_precision<Scalar>()==2) nd = 16;
      int oldprecision = str.precision(nd);

      Scalar hstep;
      if (n>1) hstep = (hmax-hmin)/(n-1);
      else hstep = ScalarFieldTraits<Scalar>::Zero();

      for (int i=0;i<n;i++)  {	
	Scalar h = hmin+i*hstep;
	x.copy(y);
	x.linComb(h,p);
	str << setprecision(6) << setw(8) << h << " "
	    << setprecision(nd) << setw(nd+6) << feval.getValue()
	    << endl;
	cerr<<h<<" "<<feval.getValue()<<endl;
      }
      str.precision(oldprecision);
    }
    catch (RVLException & e) {
      e<<"\ncalled from Scan\n";
      throw e;
    }
  }

}

#endif
