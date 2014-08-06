// lnsrchMT.H
// created by ADP and WWS
// last modified 06/09/04

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

#ifndef __RVL_ALG_LINE_SEARCH_MT
#define __RVL_ALG_LINE_SEARCH_MT

#include "alg.hh"
#include "functional.hh"
#include "terminator.hh"
#include "table.hh"
#include "lnsrch.hh"

namespace RVLUmin {
  using namespace RVLAlg;
  using namespace RVL;


  /** An implementation of the More' and Thuente line search algorithm
      (See More' and Thuente, "Line Search Algorithms with Guaranteed
      Sufficient Decrease", ACM TOMS, Vol. 20, No. 3, 286--307 (1994)).
      The description of this class follows closely that of the base
      class LineSearch, which should be consulted for details. */
  template<class Scalar>
  class LineSearchMT: public LineSearchAlg<Scalar> {

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


    Table ParamTable;

    // parameters found in ParamTable
  
    /**@name Input parameters */
    //@{
    /** (1e-4) Tolerance for decreasing function value (Armijo-Goldstein
	condition). */
    Scalar FcnDecreaseTol; 
    /** (9e-1) Slope at minimum must be smaller than SlopeDecreaseTol times the
	initial slope. */
    Scalar SlopeDecreaseTol; 
    /// (1e-20) Minimum step length 
    Scalar MinStep;
    /// (1e20) Maximum step length 
    Scalar MaxStep;
    /// (8) Maximum number of function evaluations to perform
    int MaxSample;
    /// (1e-40) Tolerance for divide-by-zero
    Scalar ZeroDivideTol;
    /** (0) DispFlag controls the amount of output sent to the screen during
	execution of the Search method.  The usual choices are
	0 --- no output
	1 --- a summary at the end of execution
	2 --- a summary at each iteration
	In addition, values of DispFlag greater than 2 may send increasing
	amount of detail to the screen (this is intended primarily for
	development use, i.e. for debugging). */
    int DispFlag;
    /** (6) DispPrecision gives the number of digits in numbers sent to
	the screen. */
    int DispPrecision;
    /** (sqrt(macheps)) Line search halts if the length of the interval of
	of uncertainty falls below IntervalTol. */
    Scalar IntervalTol;
    /// (4) factor to change mu during backtracking
    Scalar BracketIncrease;
    //@}
  
  
    /**@name Output parameters */
    //@{
    ///exit status (<0 = no reduction, 0 = success, >0 = reduction, but failed)
    mutable int TermCode;
    /// Maximum step size taken?  (1 = taken, 0 = not taken)
    mutable int MaxTkn;   
    //@}
  
    // local variables
    int OldDispPrec;
    int Itn;                   // Backtracking iteration counter
  
    ostream & str;

    /// Read parameters from parameter table and allocate temporary vectors.
    void initialize() {
      // Read the parameters
      int MinFlag = ParamTable.getValue("MinStep",MinStep);
      int MaxFlag = ParamTable.getValue("MaxStep",MaxStep);
      if (MinFlag && MaxFlag ||
	  (!MinFlag && !MaxFlag && MinStep>=MaxStep)) {
	MaxStep = 1.0e20;
	MinStep = 1.0e-20;
      }
      else if (MinFlag && !MaxFlag)
	MinStep = 1.0e-40*MaxStep;
      else if (MaxFlag && !MinFlag)
	MaxStep = 1.0e40*MinStep;
      if (ParamTable.getValue("MaxSample",MaxSample))
	MaxSample = 8;
      if (ParamTable.getValue("FcnDecreaseTol",FcnDecreaseTol))
	FcnDecreaseTol = 1e-4;
      if (ParamTable.getValue("SlopeDecreaseTol",SlopeDecreaseTol))
	SlopeDecreaseTol = 0.9;
      if (ParamTable.getValue("IntervalTol",IntervalTol))
	IntervalTol = sqrt(numeric_limits<Scalar>::epsilon());
      if (ParamTable.getValue("BracketIncrease",BracketIncrease))
	BracketIncrease = 4.0;
      if (ParamTable.getValue("ZeroDivideTol",ZeroDivideTol))
	ZeroDivideTol = 1.0e-40;
      if (ParamTable.getValue("DispFlag",DispFlag))
	DispFlag = 0;
      if (ParamTable.getValue("DispPrecision",DispPrecision))
	DispPrecision = 6;
      OldDispPrec=cout.precision(DispPrecision);
      Itn = 0;
    }
    
    /// Display results
    void displayResults(int info) {
      if (DispFlag) {
	if (TermCode==Succeed)
	  this->str<<"Line search succeeded"<<endl;
	else
	  switch (info)
	    {
	    case 2:
	      this->str<<"Line search halted: interval of uncertainty "
		"is below tolerance"<<endl;
	      break;
	    case 3:
	      this->str<<"LineSearch halted: maximum number of iterations "
		"reached"<<endl;
	      break;
	    case 4:
	      this->str<<"Line Search halted: minimum step length reached"
		  <<endl;
	      break;
	    case 5:
	      this->str<<"Line Search halted: maximum step length reached"
		  <<endl;
	      break;
	    case 6:
	      this->str<<"Line Search halted: failure to make progress; "
		"tolerances may be too small" <<endl;
	      break;
	    default:
	      this->str<<"Error in LineSearchMT::displayResults: "
		"this default case should not occur"<<endl;
	    }
	if (TermCode==NoReduction)
	  this->str<<"   Function value not reduced"<<endl;
      }
    }
  
    /// Update the interval of uncertainty and compute the next step.
    int cstep(Scalar & mux,Scalar & fx,Scalar & dgx,Scalar & muy,
	      Scalar & fy,Scalar & dgy,Scalar & mu,Scalar & fv,
	      Scalar & dg,int & bracket,Scalar & mumin,Scalar & mumax) {

      try {
	// Check the input for errors
	if (bracket && (mu<=min(mux,muy) || mu>=max(mux,muy)) ||
	    dgx*(mu-mux)>=0.0 || mumax<mumin) {
	  RVLException e;
	  e<<"Error: LineSearchMT::cstep\n";
	  e<<"error in input\n";
	  throw e;
	}
      
	// Determine if derivatives have opposite signs.
	Scalar sgnd=dg;
	if (dgx<0.0) sgnd=-dg;
	int infoc,bound;
	Scalar muf;
	int FPErr=0;
      
	// First case: a higher function value.
	// The minimum is bracketed.  If the cubic step is closer
	// to mux than the quadratic step, the cubic step is taken,
	// otherwise the average of the cubic and quadratic steps is taken.
      
	if (fv > fx) {
	  if (DispFlag>1)
	    this->str <<"   cstep: Case 1" <<endl;
	
	  infoc = 1;
	  bracket = 1;
	  bound = 1;
	  Scalar theta;
	  if (FPErr=ProtectedDivision<Scalar>(3.0*(fx-fv),
					      mu-mux,
					      theta,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) theta"
		  <<endl<<"   "<<3.0*(fx-fv)<<"/"<<mu-mux<<endl;
	    goto FPErrTarget;
	  }
	  theta = theta + dgx + dg;
	  Scalar s = max(max(abs(theta),abs(dgx)),abs(dg));
	  Scalar t1,t2,t3;
	  if (FPErr=ProtectedDivision<Scalar>(theta,
					      s,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) gamma (t1)"
		  <<endl<<"   "<<theta<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dgx,
					      s,
					      t2,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) gamma (t2)"
		  <<endl<<"   "<<dgx<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dg,
					      s,
					      t3,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) gamma (t3)"
		  <<endl<<"   "<<dg<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  Scalar gamma = s*sqrt(t1*t1-t2*t3);
	  if (mu < mux) gamma=-gamma;
	  Scalar p = (gamma-dgx)+theta;
	  Scalar q = ((gamma-dgx)+gamma)+dg;
	  Scalar r;
	  if (FPErr=ProtectedDivision<Scalar>(p,
					      q,
					      r,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) r"
		  <<endl<<"   "<<p<<"/"<<q<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muc = mux + r*(mu-mux);
	  if (FPErr=ProtectedDivision<Scalar>(fx-fv,
					      mu-mux,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) muq (t1)"
		  <<endl<<"   "<<fx-fv<<"/"<<mu-mux<<endl;
	    goto FPErrTarget;
	  }
	  t1+=dgx;
	  if (FPErr=ProtectedDivision<Scalar>(dgx,
					      t1,
					      t2,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 1) muq (t2)"
		  <<endl<<"   "<<dgx<<"/"<<t1<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muq = mux + (t2/2.0)*(mu-mux);
	  if (abs(muc-mux) < abs(muq-mux))
	    muf = muc;
	  else
	    muf = muc + (muq-muc)/2.0;
	}

	// The second case: a lower function value and derivatives of opposite
	// sign.  The minimum is bracketed.  If the cubic step is
	// closer to mux than the quadratic (secant) step, then the
	// cubic step is taken.  Otherwise the quadratic step is taken.
      
	else if (sgnd < 0.0) {
	  if (DispFlag>1)
	    this->str<<"   cstep: Case 2"<<endl;
	
	  infoc = 2;
	  bracket = 1;
	  bound = 0;
	  Scalar theta;
	  if (FPErr=ProtectedDivision<Scalar>(3.0*(fx-fv),
					      mu-mux,
					      theta,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) theta"
		  <<endl<<"   "<<3.0*(fx-fv)<<"/"<<mu-mux<<endl;
	    goto FPErrTarget;
	  }
	  theta = theta + dgx + dg;
	  Scalar s = max(max(abs(theta),abs(dgx)),abs(dg));
	  Scalar t1,t2,t3;
	  if (FPErr=ProtectedDivision<Scalar>(theta,
					      s,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) gamma (t1)"
		  <<endl<<"   "<<theta<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dgx,
					      s,
					      t2,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) gamma (t2)"
		  <<endl<<"   "<<dgx<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dg,
					      s,
					      t3,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) gamma (t3)"
		  <<endl<<"   "<<dg<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  Scalar gamma = s*sqrt(t1*t1-t2*t3);
	  if (mu > mux) gamma = -gamma;
	  Scalar p = (gamma-dg)+theta;
	  Scalar q = ((gamma-dg)+gamma)+dgx;
	  Scalar r;
	  if (FPErr=ProtectedDivision<Scalar>(p,
					      q,
					      r,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) r"
		  <<endl<<"   "<<p<<"/"<<q<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muc = mu + r*(mux-mu);
	  if (FPErr=ProtectedDivision<Scalar>(dg,
					      dg-dgx,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 2) muq (t1)"
		  <<endl<<"   "<<dg<<"/"<<dg-dgx<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muq = mu + t1*(mux-mu);
	  if (abs(muc-mu) > abs(muq-mu))
	    muf = muc;
	  else
	    muf = muq;
	}

	// Third case: a lower function value, derivatives of the same sign,
	// and the magnitude of the derivative decreases.  The cubic step
	// is only used if the cubic tends to infinity in the direction
	// of the step or if the minimum of the cubic is beyond mu.
	// Otherwise the cubic step is defined to be either mumin or mumax.
	// The quadratic (secant) step is also computed and if the minimum
	// is bracketed, then the step closest to mu x is taken, otherwise
	// the step farthest away is taken.
      
	else if (abs(dg) < abs(dgx)) {
	  if (DispFlag>1)
	    this->str<<"   cstep: Case 3"<<endl;
	
	  infoc = 3;
	  bound = 1;
	  Scalar theta;
	  if (FPErr=ProtectedDivision<Scalar>(3.0*(fx-fv),
					      mu-mux,
					      theta,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) theta"
		  <<endl<<"   "<<3.0*(fx-fv)<<"/"<<mu-mux<<endl;
	    goto FPErrTarget;
	  }
	  theta = theta + dgx + dg;
	  Scalar s = max(max(abs(theta),abs(dgx)),abs(dg));
	  Scalar t1,t2,t3;
	  if (FPErr=ProtectedDivision<Scalar>(theta,
					      s,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) gamma (t1)"
		  <<endl<<"   "<<theta<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dgx,
					      s,
					      t2,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) gamma (t2)"
		  <<endl<<"   "<<dgx<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  if (FPErr=ProtectedDivision<Scalar>(dg,
					      s,
					      t3,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) gamma (t3)"
		  <<endl<<"   "<<dg<<"/"<<s<<endl;
	    goto FPErrTarget;
	  }
	  Scalar gamma = s*sqrt(max((Scalar)0.0,(Scalar)(t1*t1-t2*t3)));
	  if (mu > mux) gamma=-gamma;
	  Scalar p = (gamma-dg)+theta;
	  Scalar q = (gamma+(dgx-dg))+gamma;
	  Scalar r;
	  if (FPErr=ProtectedDivision<Scalar>(p,
					      q,
					      r,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) r"
		  <<endl<<"   "<<p<<"/"<<q<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muc;
	  if (r<0.0 && gamma!=0.0)
	    muc = mu + r*(mux-mu);
	  else if (mu > mux)
	    muc = mumax;
	  else
	    muc = mumin;
	  if (FPErr=ProtectedDivision<Scalar>(dg,
					      dg-dgx,
					      t1,
					      ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"Protected division failed: cstep (case 3) muq (t1)"
		  <<endl<<"   "<<dg<<"/"<<dg-dgx<<endl;
	    goto FPErrTarget;
	  }
	  Scalar muq = mu + t1*(mux-mu);
	  if (bracket) {
	    if (abs(muc-mu) < abs(muq-mu))
	      muf = muc;
	    else
	      muf = muq;
	  }
	  else if (abs(mu-muc) > abs(mu-muq))
	    muf = muc;
	  else
	    muf = muq;
	}
      
	// Fourth case. a lower function value, derivatives of the same
	// sign, and the magnitude of the derivative does not decrease.
	// If the minimum is not bracketed, the step is either mumin or
	// mumax.  Otherwise the cubic step is taken.
      
	else {
	  if (DispFlag>1)
	    this->str<<"   cstep: Case 4"<<endl;
	
	  infoc = 4;
	  bound = 0;
	  if (bracket) {
	    Scalar theta;
	    if (FPErr=ProtectedDivision<Scalar>(3.0*(fv-fy),
						muy-mu,
						theta,
						ZeroDivideTol)) {
	      if (DispFlag>1)
		this->str<<"Protected division failed: cstep (case 4) theta"
		    <<endl<<"   "<<3.0*(fv-fy)<<"/"<<muy-mu<<endl;
	      goto FPErrTarget;
	    }
	    theta = theta + dgy + dg;
	    Scalar s = max(max(abs(theta),abs(dgy)),abs(dg));
	    Scalar t1,t2,t3;
	    if (FPErr=ProtectedDivision<Scalar>(theta,
						s,
						t1,
						ZeroDivideTol)) {
	      if (DispFlag>1)
		this->str<<"Protected division failed: cstep (case 4) gamma (t1)"
		    <<endl<<"   "<<theta<<"/"<<s<<endl;
	      goto FPErrTarget;
	    }
	    if (FPErr=ProtectedDivision<Scalar>(dgy,
						s,
						t2,
						ZeroDivideTol)) {
	      if (DispFlag>1)
		this->str<<"Protected division failed: cstep (case 4) gamma (t2)"
		    <<endl<<"   "<<dgy<<"/"<<s<<endl;
	      goto FPErrTarget;
	    }
	    if (FPErr=ProtectedDivision<Scalar>(dg,
						s,
						t3,
						ZeroDivideTol)) {
	      if (DispFlag>1)
		this->str<<"Protected division failed: cstep (case 4) gamma (t3)"
		    <<endl<<"   "<<dg<<"/"<<s<<endl;
	      goto FPErrTarget;
	    }
	    Scalar gamma = s*sqrt(t1*t1-t2*t3);
	    if (mu > muy) gamma=-gamma;
	    Scalar p = (gamma-dg)+theta;
	    Scalar q = ((gamma-dg)+gamma)+dgy;
	    Scalar r;
	    if (FPErr=ProtectedDivision<Scalar>(p,
						q,
						r,
						ZeroDivideTol)) {
	      if (DispFlag>1)
		this->str<<"Protected division failed: cstep (case 4) r"
		    <<endl<<"   "<<p<<"/"<<q<<endl;
	      goto FPErrTarget;
	    }
	    Scalar muc = mu + r*(muy-mu);
	    muf = muc;
	  }
	  else if (mu > mux)
	    muf = mumax;
	  else
	    muf = mumin;
	}
      
      FPErrTarget: ;
      
	// Update the interval of uncertainty.  This update does not
	// depend on the new step or the case analysis above.
      
	if (fv > fx) {
	  muy = mu;
	  fy = fv;
	  dgy = dg;
	}
	else {
	  if (sgnd < 0.0) {
	    muy = mux;
	    fy = fx;
	    dgy = dgx;
	  }
	  mux = mu;
	  fx = fv;
	  dgx = dg;
	}

	// Bisect the interval if there was a floating point error in
	// the calculation of muf
      
	if (FPErr && bracket)
	  muf = 0.5*(mux+muy);
	else if (FPErr)
	  muf = mumax;
      
	// Compute the new step and safeguard it.
      
	muf = min((Scalar)mumax,(Scalar)muf);
	muf = max((Scalar)mumin,(Scalar)muf);
	mu = muf;
	if (bracket && bound)
	  if (muy > mux)
	    mu = min((Scalar)(mux+0.66*(muy-mux)),(Scalar)mu);
	  else
	    mu = max((Scalar)(mux+0.66*(muy-mux)),(Scalar)mu);
	return infoc;
      }
      catch(RVLException & e) {
	e<<"\ncalled from LineSearchMT::cstep\n";
	throw e;
      }
    }
  
    // copy constructors --- disabled
    LineSearchMT(const LineSearchMT<Scalar> &);
  
  public:

    /** Initialize a line search with the info for the 
	parameter table.
    */
    LineSearchMT( const Space<Scalar> & sp, 
		  string parname = "",
		  string prefix = "LineSearchMT",
		  ostream & _str = cout) 
      : LineSearchAlg<Scalar>(sp), ParamTable(parname,prefix), str(_str) {
      try { initialize(); }
      catch (RVLException & e) {
	e<<"\ncalled from LineSearchMT constructor\n";
	throw e;
      }
    }

    // Destructor
    virtual ~LineSearchMT() {}

    bool query() {if (TermCode) return false; return true; }
    
    void run() {
      try {
      
	if (DispFlag)
	  this->str<<"===============LineSearchMT==============="<<endl;
      
	// Compute the length of the Newton step, and scale
	// the search direction (if necessary)
      
	const Vector<Scalar> & x0 = this->getBasePoint();
	Vector<Scalar> & dx = this->getSearchDirection();
	
	FunctionalEvaluation<Scalar> & feval = this->getFunctionalEvaluation();
	Vector<Scalar> & x = feval.getPoint();

	Scalar newtlen=this->getBaseGradient().norm();
      
	Scalar t = feval.getMaxStep(dx);
	Scalar fmaxs;
	if (t <= 0.0) {
	  if (DispFlag)
	    this->str<<"Error in LineSearchMT: vector not in domain"
	      " (fmaxs == 0.0)"<<endl;
	  
	  TermCode = NoReduction;
	  displayResults(6);
	  return;
	}
	if (t==numeric_limits<Scalar>::max())
	  fmaxs = t;
	else
	  fmaxs = 0.99*newtlen*t;
	Scalar fcnmaxstep = min(fmaxs,MaxStep);
      
	if (newtlen > fcnmaxstep) {   // ||dir|| > MaxStep
	  Scalar scale;
	  if (ProtectedDivision(fcnmaxstep,
				newtlen,
				scale,
				ZeroDivideTol)) {
	    if (DispFlag>1)
	      this->str<<"ProtectedDivision error: newtlen too small" 
		  <<"   "<<fcnmaxstep<<"/"<<newtlen<<endl;
	    TermCode = NoReduction;
	    displayResults(6);
	    return;
	  }
	  // scale descent direction
	  dx.scale(scale);
	  newtlen = fcnmaxstep;
	  if (DispFlag>1)
	    this->str<<"Vector dir is too long; multiply by "<<scale<<endl;
	}
      
	// Compute the initial slope; exit if not negative
         
	Scalar dginit = this->getBaseGradient().inner(dx);
	Scalar fcur = feval.getValue();
	if (dginit >= 0.0) {
	  if (DispFlag) {
	    this->str<<"Initial slope nonnegative in line search"<<endl;
	    this->str<<"dginit: "<<dginit<<endl;
	  }
	  TermCode = InitialSlopePositive;
	  displayResults(6);
	  return;
	}
	Scalar MinMu,MaxMu;
	if (ProtectedDivision<Scalar>(MinStep,
				      newtlen,
				      MinMu,
				      ZeroDivideTol) ||
	    ProtectedDivision<Scalar>(fcnmaxstep,
				      newtlen,
				      MaxMu,
				      ZeroDivideTol)) {
	  if (DispFlag>1)
	    this->str<<"Protected division failed (computing MinMu and MaxMu)"
		<<"   "<<MinStep<<"/"<<newtlen<<", "
		<<fcnmaxstep<<"/"<<newtlen<<endl;
	  TermCode = NoReduction;
	  displayResults(6);
	  return;
	}
      
	if (DispFlag>1) {
	  Scalar xnorm=x.norm();
	  this->str<<"f = "<<fcur<<endl;
	  this->str<<"||x|| = "<<xnorm<<endl;
	  this->str<<"||dir|| = "<<newtlen<<endl;
	  this->str<<"Initial slope = "<<dginit<<endl;
	  this->str<<"Maximum step to boundary: "<<fmaxs<<endl;
	  this->str<<"Minimum mu: "<<MinMu<<endl;
	  this->str<<"Maximum mu: "<<MaxMu<<endl;
	}

	int bracket = 0;
	int stage1 = 1;
	int NumSamples = 0;
	Scalar finit = fcur;
	Scalar dgtest = FcnDecreaseTol*dginit;
	Scalar width = MaxStep - MinStep;
	Scalar width1 = width*2.0;
      
	Scalar mux = 0.0;
	Scalar fx = finit;
	Scalar dgx = dginit;
	Scalar muy = 0.0;
	Scalar fy = finit;
	Scalar dgy = dginit;
	Scalar mu = 1.0;
	mu = min(mu,MaxMu);
	mu = max(mu,MinMu);
	Scalar mumin,mumax;
      
	MaxTkn = 0;
	int infoc;
	Scalar dg;
      
	Scalar fnext;

	while(true) {
	  // Set the minimum and maximum steps to correspond to the present
	  // interval of uncertainty.
	
	  if (bracket) {
	    mumin = min(mux,muy);
	    mumax = max(mux,muy);
	  }
	  else {
	    mumin = mux;
	    mumax = mu + BracketIncrease*(mu-mux);
	  }
	  
	  // Force the step to be within the bounds
	  mu = max(mu,MinMu);
	  mu = min(mu,MaxMu);

	  // In case of failure, let mu correspond to the best point
	  // Note: This is slightly different from the More'-Thuente
	  // code.  In case the interval of uncertainty becomes too
	  // small, we accept the new mu as the "best point", since
	  // any point in the small interval should equally acceptable,
	  // and it it likely that the new mu is better than the
	  // previous best point.
	
	  if ((bracket && (mu <= mumin || mu >= mumax)) ||
	      NumSamples >= MaxSample-1)
	    mu = mux;
	
	  // Evaluate the function and gradient at mu - enclose in block
	  // to localize fevalnext, which is not needed apart from this
	  x.copy(x0);
	  x.linComb(mu,dx);
	
	  fnext = feval.getValue();  
	  NumSamples++;
	  dg = feval.getGradient().inner(dx);
	  Scalar ftest1 = finit + mu*dgtest;
	  if (DispFlag>1)
	    this->str<<NumSamples<<" mu = "<<mu<<" f = "<<fnext<<" slope = " 
		<<dg<<endl;
	
	  // Check for convergence
	
	  int info=0;
	  if ((bracket && (mu <= mumin || mu >= mumax)))
	    info = 6;
	  else if (mu == MaxMu && fnext <= ftest1 && dg <= dgtest) {
	    info = 5;
	    MaxTkn = 1;
	  }
	  else if (mu == MinMu && (fnext > ftest1 || dg >= dgtest))
	    info = 4;
	  else if (NumSamples >= MaxSample)
	    info = 3;
	  else if (bracket && mumax-mumin <= IntervalTol*mumax)
	    info = 2;
	
	  if (info) {
	    if (fnext >= finit)
	      TermCode = NoReduction;
	    else if (info == 3)
	      TermCode = TooManyFunctionEvals;
	    else
	      TermCode = Fail;    
	    displayResults(info);
	    return;
	  }
	  if (fnext <= ftest1 && abs(dg) <= SlopeDecreaseTol*(-dginit)) {
	    TermCode = Succeed;
	    displayResults(1);
	    
	    this->step = mu;
	    return;
	  }
	
	  // In the first stage, we seek a step for which the modified
	  // function has a nonpositive value and nonnegative derivative.
	
	  if (stage1 && fnext <= ftest1 &&
	      dg >= min(FcnDecreaseTol,SlopeDecreaseTol)*dginit)
	    stage1 = 0;
	
	  // We use the modified function to predict the step only if
	  // we do not have a step for which the modified function has
	  // a nonpositive function value and nonnegative derivative,
	  // and if a lower function value has been obtained but the
	  // decrease is not sufficient.
	
	  if (stage1 && fnext <= fx && fnext >= ftest1) {

	    // Define the modified function and derivative values.
	    Scalar fm = fnext - mu*dgtest;
	    Scalar fxm = fx - mux*dgtest;
	    Scalar fym = fy - muy*dgtest;
	    Scalar dgm = dg - dgtest;
	    Scalar dgxm = dgx - dgtest;
	    Scalar dgym = dgy - dgtest;
	  
	    // Update the interval of uncertainty and compute the next step.  
	    infoc = cstep(mux,fxm,dgxm,muy,fym,dgym,mu,fm,dgm,bracket,
			  mumin,mumax);
	    
	    // Reset the function and gradient values.
	    fx = fxm + mux*dgtest;
	    fy = fym + muy*dgtest;
	    dgx = dgxm + dgtest;
	    dgy = dgym + dgtest;
	  }
	  else
	    // Update the interval of uncertainty and compute the new step.
	    infoc = cstep(mux,fx,dgx,muy,fy,dgy,mu,fnext,dg,bracket,mumin,mumax);
	
	  // Force a sufficient decrease in the size of the interval of
	  // uncertainty.
	
	  if (bracket) {
	    if (abs(muy-mux) >= 0.66*width1)
	      mu = mux + 0.5*(muy-mux);
	    width1 = width;
	    width = abs(muy-mux);
	  }
	}
      }
      catch(RVLException & e) {
	e<<"\ncalled from LineSearchMT::search\n";
	throw e;
      }
    }

    // Access to parameter table. 
    Table & accessParameters() { return ParamTable; }  

    /// write methods to print out useful information about the object.
    void write(RVLException & e) {
      e<<"LineSearchMT More-Thuente line search algorithm\n";
      //e<<ParamTable.write(e);
    }
    ostream & write(ostream & str) {
      str<<"LineSearchMT More-Thuente line search algorithm\n";
      ParamTable.write(str);
      return str;
    }
  };

}

#endif  
