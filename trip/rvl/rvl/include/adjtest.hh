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

#ifndef __RVL_ADJTEST
#define __RVL_ADJTEST

#include "op.hh"

namespace RVL {

  /** Test of adjoint relationship between applyOp and applyAdjOp
      methods of a LinearOp. constructs random vectors x in domain,
      y in range, applies image method to x, adjoint image method to
      y, and compares inner products. */
  
  template<typename Scalar>
  bool AdjointTest(LinearOp<Scalar> const & op,
		   FunctionObject & randomize,
		   ostream & str,
		   int tol=100) {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

    try {

      Vector<Scalar> xin(op.getDomain());
      Vector<Scalar> yin(op.getRange());
      Vector<Scalar> xout(op.getDomain());
      Vector<Scalar> yout(op.getRange());

      xin.eval(randomize);
      yin.eval(randomize);
  
      xout.zero();
      yout.zero();
  
      op.applyOp(xin,yout);
      op.applyAdjOp(yin,xout);
      Scalar ipdom=xin.inner(xout);
      Scalar iprng=yout.inner(yin);
      atype denom=yout.norm()*yin.norm();
      atype eps = tol*numeric_limits<atype>::epsilon();
      // here |Ax| is used as an estimate for |A||x|
	
      int nd;
      if (numeric_precision<atype>()==1) nd = 8;
      if (numeric_precision<atype>()==2) nd = 16;
      str<<setprecision(nd)<<setiosflags(ios_base::scientific);

      atype aip=abs(ipdom-iprng);
      if (aip<eps*denom) {
	str<<"AdjointTest"<<endl<<endl;
	str<<"adjoint relation holds:";
	atype rat;
	if (ProtectedDivision<atype>(aip,denom,rat)) {
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<Ax,y>  = "<<iprng<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<x,A'y> = "<<ipdom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|Ax||y|  = "<<denom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|<Ax,y>-<x,A'y>| = "<<aip<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"< 100*macheps*|Ax||y| = "<<eps*denom<<endl;
	}
	else {
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|<Ax,y>-<x,A'y>|/|Ax||y| = "<<rat<<endl;
	  str<<setprecision(nd)<<setiosflags(ios::scientific)<<"< 100*macheps = "<<eps<<endl;
	  str<<"<Ax,y>  = "<<setprecision(nd)<<setiosflags(ios::scientific)<<iprng<<endl;
	  str<<setprecision(nd)<<setiosflags(ios::scientific)<<"<x,A'y> = "<<ipdom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios::scientific)<<"|Ax||y|  = "<<denom<<endl;
	}
	return true;
      }
      else {
	str<<"AdjointTest"<<endl;
	str<<"adjoint relation fails: ";
	
	atype rat;
	if (ProtectedDivision<atype>(aip,denom,rat)) {
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<Ax,y>  = "<<iprng<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<x,A'y> = "<<ipdom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|Ax||y|  = "<<denom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|<Ax,y>-<x,A'y>| = "<<aip<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<">= 100*macheps*|Ax||y| = "<<eps*denom<<endl;     
	}
	else {
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|<Ax,y>-<x,A'y>|/|Ax||y| = "<<rat<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<">= 100*macheps = "<<eps<<endl;     
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<Ax,y>  = "<<iprng<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"<x,A'y> = "<<ipdom<<endl;
	  str<<setprecision(nd)<<setiosflags(ios_base::scientific)<<"|Ax||y|  = "<<denom<<endl;
	}
      }
      return false;
    }
    catch (RVLException & e) {
      e<<"\ncalled from AdjointTest\n";
      throw e;
    }
  }
}

#endif
