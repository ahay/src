// lnsrchDS.H
// created by ADP and WWS
// last modified 08.12.04

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

#ifndef __RVL_ALG_LINE_SEARCH_DS
#define __RVL_ALG_LINE_SEARCH_DS

#include "umintable.hh"
#include "lnsrch.hh"

namespace RVLUmin {
  using namespace RVLAlg;
  using namespace RVL;

  /** An implementation of a backtracking line search algorithm using 
      cubic interpolation (See Dennis and Schnabel, "Numerical Methods
      for Unconstrained Optimization and Nonlinear Equations",
      Prentice-Hall (1983)). 
      The description of this class differs from the description of the
      base class LineSearch in only one regard: the second
      criteria for defining an acceptable step is 
      \f$
      \phi'(a) \ge tol \phi'(0).
      \f$
      This prevents short steps by requiring that the slope be increased
      from its initial (negative) value, although it does not require
      that the step reach an approximate minimizer.  This criterion is
      sufficient to allow a positive definite update in a BFGS algorithm.

      Adapted from Dennis-Schnabel line search implementation in HCL1.0,
      authored by Mark Gockenbach.
    
      ADP and WWS, Spring 04.
  */
  template<class Scalar>
  class LineSearchDS: public LineSearchAlg<Scalar> {

  private:

    /**@name Term codes */
    //@{
    enum {
      NoErr = 0,
      Succeed = 0,
      NoReduction = -2,
      InitialSlopePositive = -1,
      TooManyFunctionEvals = 2,
      Fail = 3
    } ErrorCodes;
    //@}

    // parameters found in UMinTable - see docs on that class for
    // meaning and defaults 

    Scalar FcnDecreaseTol; 
    Scalar SlopeDecreaseTol;
    int MaxSample;
    Scalar MinStep;
    Scalar MaxStep;
    Scalar ZeroDivideTol;
    int DispFlag;
    int DispPrecision;
    Scalar BracketIncrease;

    /**@name Output parameters */
    //@{
    ///exit status  (<0 = no reduction, 0 = success, >0 = reduction, but failed)
    mutable int TermCode;
    ///maximum step size taken?  (1 = taken, 0 = not taken)
    mutable int MaxTkn;   
    //@}
   
    /// Read parameters from parameter table
    void initialize(UMinTable<Scalar> & holder) {

      if (holder.getValue("MinStep",MinStep) ||
	  holder.getValue("MaxStep",MaxStep) ||
	  holder.getValue("MaxSample",MaxSample) ||
	  holder.getValue("FcnDecreaseTol",FcnDecreaseTol) ||
          holder.getValue("SlopeDecreaseTol",SlopeDecreaseTol) ||
	  holder.getValue("BracketIncrease",BracketIncrease) ||
	  holder.getValue("ZeroDivideTol",ZeroDivideTol) ||
	  holder.getValue("DispFlag",DispFlag) ||
	  holder.getValue("DispPrecision",DispPrecision)) {
	RVLException e;
	e<<"Error: LineSearchDS::initialize\n";
	e<<"UMinTable object corrupted\n";
	holder.write(cerr);
	throw e;
      }
    }

    ostream & str;

    LineSearchDS();

  public:

    /** Initialize internal data from a table loaded from a file.
    */
    LineSearchDS(const Space<Scalar> & sp, 
		 UMinTable<Scalar> & holder,
		 ostream & _str = cout) 
      : LineSearchAlg<Scalar>(sp), str(_str) {
      try { initialize(holder); }
      catch (RVLException & e) {
	e<<"\ncalled from LineSearchDS constructor\n";
	throw e;
      }
    }

    /** copy constructor */
    LineSearchDS(const LineSearchDS<Scalar> & a) 
      : LineSearchAlg<Scalar>(a),
	DispFlag(a.DispFlag),
	MaxSample(a.MaxSample),
	FcnDecreaseTol(a.FcnDecreaseTol),
	SlopeDecreaseTol(a.SlopeDecreaseTol),
	BracketIncrease(a.BracketIncrease),
	MinStep(a.MinStep),
	MaxStep(a.MaxStep),
	ZeroDivideTol(a.ZeroDivideTol),
	DispPrecision(a.DispPrecision),
	str(a.str)
    {}

    // Destructor
    virtual ~LineSearchDS() {}

    bool query() {if (TermCode) return false; return true; }
  
    void run() {
      try {
	int NumFcnSampled = 0;
	TermCode = NoErr;
	MaxTkn = 0;

	if (DispFlag)
	  str<<"===============LineSearchDS==============="<<endl;

	// Compute the length of the Newton step, and scale
	// the search direction (if necessary)
	const Vector<Scalar> & x0 = this->getBasePoint();
	Vector<Scalar> & x = this->getTrialPoint();
	Vector<Scalar> & dx = this->getSearchDirection();
	FunctionalEvaluation<Scalar> & fx = LineSearchAlg<Scalar>::getFunctionalEvaluation();

	// no need to copy!!!!
	//	x.copy(x0);
	Scalar newtlen=dx.norm();
      
	Scalar t = fx.getMaxStep(dx);
	if (t<=0.0) {
	  RVLException e;
	  e<<"Error: LineSearchDS::run\n";
	  e<<"any positive step in direction indicated is infeasible\n";
	  if (DispFlag) e.write(str);
	  throw e;
	}
      
	Scalar fmaxs;
	if (t == numeric_limits<Scalar>::max())
	  fmaxs = t;
	else
	  fmaxs = 0.99*newtlen*t;
      
	Scalar fcnmaxstep = min(fmaxs,MaxStep);

	if (newtlen>fcnmaxstep) {   // ||dir|| > MaxStep
	  Scalar scale;
	  if (ProtectedDivision<Scalar>(fcnmaxstep,
					newtlen,
					scale,
					ZeroDivideTol)) {
	    RVLException e;
	    e<<"Error: LineSearchDS::run\n";	    
	    e<<"ProtectedDivision error: newtlen too small";
	    <<"   "<<fcnmaxstep<<"/"<<newtlen<<endl;
	    if (DispFlag) e.write(str);
	    throw e;
	  }

	  // scale descent direction
	  dx.scale(scale);
	  newtlen = fcnmaxstep;
	  if (DispFlag > 1)
	    str<<"Vector dir is too long; multiply by "
		<<scale<<endl;
	}
      
	// Compute the initial slope; exit if not negative
      
	Scalar initslope;
	Scalar fcur, fnext;
      
	fcur = fx.getValue();
	//	cerr<<"fx.getGradient in LineSearchDS"<<endl;
	initslope = fx.getGradient().inner(dx);

	if (initslope>=0.0) {
	  if (DispFlag) {
	    str<<"Initial slope nonnegative in line search"<<endl;
	    str<<"initslope: "<<initslope<<endl;
	  }
	  TermCode=InitialSlopePositive;
	  return;
	}
      
	Scalar minmu,maxmu;
	if (ProtectedDivision<Scalar>(MinStep,
				      newtlen,
				      minmu,
				      ZeroDivideTol)) 
	  minmu=numeric_limits<Scalar>::min();

	if (ProtectedDivision<Scalar>(fcnmaxstep,
				      newtlen,
				      maxmu,
				      ZeroDivideTol)) 
	  maxmu=numeric_limits<Scalar>::max();

	if (DispFlag > 1) {
	  Scalar xnorm=x.norm();
	  str<<"f(x) = "<<fcur<<endl;
	  str<<"||x|| = "<<xnorm<<endl;
	  str<<"||dir|| = "<<newtlen<<endl;
	  str<<"Initial slope = "<<initslope<<endl;
	  str<<"Maximum step to boundary: "<<fmaxs<<endl;
	  str<<"Minimum mu: "<<minmu<<endl;
	  str<<"Maximum mu: "<<maxmu<<endl;
	}
      
	// Initialize mu and enter the main loop
      
	Scalar mu=1.0;
	mu = min(mu,maxmu);
	mu = max(mu,minmu);
	Scalar mu0=mu;

	do {
	  // Declare some variables used in the iteration
	
	  Scalar muprev,mudiff,muincr,mutemp,mulo;
	  Scalar fhi,flo,fprev;
	  Scalar newslope;
	  Scalar t1,t2,t3,a,b,disc;
	  Scalar v1,v2,w1,w2;       // eric
	
	  // get next value and gradient

	  if (DispFlag > 1)
	    str<<"mu = "<<mu<<endl;
	  //	  cerr<<"update x\n";
	  x.copy(x0);
	  x.linComb(mu,dx);
	  fnext = fx.getValue();
	  NumFcnSampled++;

	  // Check for sufficient decrease in the function value
	
	  if (fnext <= fcur + FcnDecreaseTol*mu*initslope) {  
	  
	    // xnext satisfactory

	    if (DispFlag > 1)
	      str<<NumFcnSampled<<" Sufficient decrease in f; f+ = "
		  <<fnext<<endl;
	  
	    // The function value has decreased sufficiently; compute
	    // the gradient and check to see if the
	    // slope in the search direction has increased
	    // sufficiently (this prevents very short steps)
	  
	    //	    cerr<<"fx.getGradient in LineSearchDS"<<endl;
	    newslope = fx.getGradient().inner(dx);
	  
	    if (newslope < SlopeDecreaseTol*initslope) {
	    
	      // The slope in the search direction has not increased
	      // sufficiently.

	      if (DispFlag > 1)
		str<<"Insufficient decrease in slope; new slope = "
		    <<newslope<<endl;
	    
	      if (mu==mu0 && mu<maxmu) {
	      
		// This is the first step.  Keep increasing the step 
		// length until there is simultaneously sufficient
		// decrease in the function and sufficient increase
		// in the gradient, or until the maximum step length
		// is reached or until there is no longer sufficient
		// decrease in the function.
	      
		do {
		  muprev = mu;
		  fprev = fnext;
		  if (BracketIncrease*mu >= maxmu)
		    mu=0.99*maxmu;
		  else
		    mu *= BracketIncrease;
		  //		  cerr<<"update x"<<endl;
		  x.copy(x0);
		  x.linComb(mu,dx);
		  fnext = fx.getValue();
		  NumFcnSampled++;

		  if (DispFlag > 1)
		    str<<NumFcnSampled<<" Increase mu to "<<mu 
			<<" f+ = "<<fnext<<endl;

		  if (fnext <= fcur + FcnDecreaseTol*mu*initslope) {
		    //		    cerr<<"fx.getGradient in LineSearchDS\n";
		    newslope = fx.getGradient().inner(dx);
		    if (DispFlag > 1)
		      str<<"Sufficient decrease in function value; "
			"new slope = "<<newslope<<endl;
		  }
		} while((fnext <= fcur +FcnDecreaseTol*mu*initslope &&
			 newslope < SlopeDecreaseTol*initslope &&
			 mu < 0.99*maxmu && NumFcnSampled < MaxSample)
			|| fnext == FLT_MAX);
	      }
	    
	      if (mu<mu0 || (mu>mu0 &&
			     fnext>fcur+FcnDecreaseTol*mu*initslope)) {
	      
		// There is an acceptable step length between mulo and the
		// muhi; find it by using successive quadratic interpolation.
	      
		mulo = min(mu,muprev);
		mudiff = abs(mu-muprev);
		if (mu<mulo) {
		  flo=fnext;
		  fhi=fprev;
		}
		else {
		  flo=fprev;
		  fhi=fnext;
		}

		if (DispFlag > 1)
		  str<<"An acceptable step length has been bracketed: ["
		      <<mulo<<","<<mulo+mudiff<<"]"<<endl;

		do {
		  if (ProtectedDivision<Scalar>(-newslope*mudiff*mudiff,
						2.0*(fhi-(flo+newslope*mudiff)),
						muincr,
						ZeroDivideTol)) {
		    if (DispFlag > 1)
		      str<<"Protected division failed: searching "
			  <<"within bracket"<<endl
			  <<"   "<<-newslope*mudiff*mudiff<<"/"
			  <<2.0*(fhi-(flo+newslope*mudiff))<<endl;
		    muincr = 0.2*mudiff;
		  }       
		  if (muincr<0.2*mudiff)
		    muincr=0.2*mudiff;

		  mu = mulo + muincr;
		  //		  cerr<<"update x"<<endl;
		  x.copy(x0);
		  x.linComb(mu,dx);
		  fnext = fx.getValue();
		  NumFcnSampled++;

		  if (DispFlag > 1)
		    str<<NumFcnSampled 
			<<" Quadratic interpolation: mu = "<<mu 
			<<" f+ = "<<fnext<<endl;
		
		  if (fnext > fcur+FcnDecreaseTol*mu*initslope) {
		    // insufficient decrease in function value
		    mudiff=muincr;
		    fhi=fnext;
		    if (DispFlag > 1)
		      str<<" Insufficient decrease in function value" 
			  <<endl;
		  }
		  else {
		    if (DispFlag > 1)
		      str<<"Sufficient decrease in function value "
			  <<endl;
		    //		    cerr<<"fx.getGradient in LineSearchDS\n";
		    newslope = fx.getGradient().inner(dx);
		    if (DispFlag > 1)
		      str<<"new slope = "<<newslope<<endl;

		    if (newslope < SlopeDecreaseTol*initslope) {
		      mulo=mu;
		      mudiff=mudiff-muincr;
		      flo=fnext;
		      if (DispFlag > 1)
			str<<"Insufficient decrease in slope; "
			  "mudiff = "<<mudiff<<endl;
		    }
		  }
		} while(newslope<SlopeDecreaseTol*initslope 
			&& mudiff>=minmu && 
			NumFcnSampled<MaxSample);
	      }
	      if (newslope < SlopeDecreaseTol*initslope) {
		TermCode=Fail;
		if (DispFlag) {
		  str<<NumFcnSampled<<" mu = "<<mu<<endl;
		  str<<"Failed to satisfy SlopeDecreaseTol condition"
		      <<endl;
		}
		if (mu*newtlen > 0.99*MaxStep)
		  MaxTkn = 1;
		
		this->step = mulo;
		
		return;
	      }
	      else {
		if (DispFlag > 1)
		  str<<"Sufficient decrease in gradient; newslope = "
		      <<newslope<<endl;
	      }
	    }
	    else {
	      if (DispFlag > 1)
		str<<"Sufficient decrease in gradient; newslope = " 
		    <<newslope<<endl;
	    }
	    if (mu*newtlen > 0.99*MaxStep)
	      MaxTkn = 1;
	  
	    // The line search was successful, so return xnext
	    if (DispFlag)
	      str<<"Number of fcn evals = "<<NumFcnSampled<<" rel step = "<<mu<<endl;
	    
	    TermCode=Succeed;
	    this->step = mu;
	    return;
	  } // end (huge) if block
	  else {
	    // failed to decrease function value sufficiently
	    if (DispFlag > 1)
	      str<<NumFcnSampled
		  <<" Insufficient decrease in function value" 
		  <<" f+ = "<<fnext<<endl;
	  
	    // Check to see if mu is too small
	  
	    if (mu<=minmu) {
	    
	      if (fnext < fcur) {
		TermCode = Fail;
		if (DispFlag) {
		  str<<NumFcnSampled<<" mu = "<<mu<<endl;
		  str<<"Failed to make sufficient progress (mu too small)"
		      <<endl;
		}
		return;
	      }
	      else {
		TermCode=NoReduction;
		if (DispFlag) {
		  str<<NumFcnSampled<<" mu = "<<mu<<endl;
		  str<<"Failed to make progress (mu too small) "<<endl;
		}
		return;
	      }
	    }
	  
	    // Otherwise, perform cubic backtrack (quadratic on
	    // first step).
	  
	    else {
	      if (mu == mu0) {
		if (ProtectedDivision<Scalar>(-initslope,
					      2.0*(fnext-fcur-mu*initslope),
					      mutemp,
					      ZeroDivideTol)) {
		  if (DispFlag > 1)
		    str<<"Protected division failed: quadratic "
			<<"interpolation"<<endl
			<<"   "<<-newslope*mudiff*mudiff<<"/"
			<<2.0*(fhi-(flo+newslope*mudiff))<<endl;
		  mutemp = 0.5*mu;
		}
		if (DispFlag > 1)
		  str<<"Quadratic backtrack: mutemp = "<<mutemp<<endl;
	      }
	      else {
		t1 = fnext-fcur-mu*initslope;
		t2 = fprev-fcur-muprev*initslope;
	      
		int FPErr = 0;
		if (FPErr=ProtectedDivision<Scalar>(1.0,
						    mu-muprev,
						    t3,
						    ZeroDivideTol)) {
		  if (DispFlag > 1)
		    str<<"Protected division failed: cubic "
		      "interpolation (t3)"<<endl
			<<1.0<<"/"<<mu-muprev<<endl;
		  goto FPErrTarget;
		}
		if (FPErr=ProtectedDivision<Scalar>(t1,
						    mu*mu,
						    v1,
						    ZeroDivideTol)) {
		  if (DispFlag > 1)
		    str<<"Protected division failed: cubic "
		      "interpolation (v1)"<<endl
			<<t1<<"/"<<mu*mu<<endl;
		  goto FPErrTarget;
		}
		if (FPErr=ProtectedDivision<Scalar>(t2,
						    muprev*muprev,
						    v2,
						    ZeroDivideTol)) {
		  if (DispFlag > 1)
		    str<<"Protected division failed: cubic "
		      "interpolation (v2)"<<endl
			<<t2<<"/"<<muprev*muprev<<endl;
		  goto FPErrTarget;
		}
		a = t3*(v1-v2);
		b = t3*(mu*v2-muprev*v1);
		disc = b*b-3.0*a*initslope;
		if( abs(a) < ZeroDivideTol ) {
		  if (FPErr=ProtectedDivision<Scalar>(-initslope,
						      2.0*b,
						      w1,
						      ZeroDivideTol)) {
		    if (DispFlag > 1)
		      str<<"Protected division failed: cubic "
			"interpolation (w1)"<<endl
			  <<-initslope<<"/"<<2.0*b<<endl;
		    goto FPErrTarget;
		  }	
		  mutemp = w1;
		}      
		else { 
		  if (FPErr=ProtectedDivision<Scalar>(sqrt(disc)-b,
						      3.0*a,
						      w2,
						      ZeroDivideTol)) {
		    if (DispFlag > 1)
		      str<<"Protected division failed: cubic "
			"interpolation (w2)"<<endl
			  <<sqrt(disc)-b<<"/"<<3.0*a<<endl;
		    goto FPErrTarget;
		  }
		  mutemp = w2;
		}
	      FPErrTarget:
		if (FPErr)
		  mutemp = 0.5*mu;
		else
		  mutemp=(a==0.0?w1:w2);
	      
		if (DispFlag > 1)
		  str<<"Cubic backtrack: mutemp = "<<mutemp<<endl;
	      }
	      mutemp = max(minmu,min((Scalar)0.5*mu,(Scalar)mutemp));
	    }
	    muprev = mu;
	    fprev = fnext;
	    mu = max((Scalar)0.1*mu,(Scalar)mutemp);
	  }
	} while (NumFcnSampled < MaxSample);
      
	if (fnext<fcur) {
	  TermCode = TooManyFunctionEvals;
	  if (DispFlag)
	    str<<NumFcnSampled<<" mu = "<<mu<<endl;
	  return;
	}
	else {
	  TermCode = NoReduction;
	  if (DispFlag)
	    str<<NumFcnSampled<<" mu = "<<mu<<endl;
	  return;
	}
      }
      catch(RVLException & e) {
	e<<"\ncalled from LineSearchDS::search\n";
	throw e;
      }
    }

    // write methods to print out useful information about the object.
    void write(RVLException & e) {
      e<<"LineSearchDS Dennis-Schnabel line search algorithm\n";
    }
    ostream & write(ostream & str) {
      str<<"LineSearchDS Dennis-Schnabel line search algorithm\n";
      return str;
    }
  };
 
}

#endif
