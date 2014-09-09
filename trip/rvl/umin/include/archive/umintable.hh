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

// author: WWS
// date of first draft: 10.12.04

#ifndef __RVLALG_UMIN_LBFGSTBL_H
#define __RVLALG_UMIN_LBFGSTBL_H

/** Table object augmented with default constructor which tests for presence
    values for keys required in constructing various update and line search 
    objects, supplies default values for any absent keys.

    Keys defined by UMinTable: [format: type name: description (default value)]
    
    string UpdateAlg: tag for update method. Allowable values:
        LBFGS = limited memory BFGS (LBFGS)

    string LineSearchAlg: tag for line search method. Allowable values:
        DS = DennisSchnabel,
        BT = simple backtracking line search  (BT)

    int MaxItn: Maximum number of iterations.(100) 

    Scalar Typf: Typical value of the functional near the solution.
	Setting	this makes the gradient stopping tolerance more
	meaningful. (1.0) 

    Scalar TypxNorm: Typical norm of the unknown vector x near the solution.
	Setting this makes the gradient stopping tolerance more
	meaningful. (1.0) 

    Scalar GradTol: Gradient tolerance. The algorithm attempts to locate
	a point where the relative gradient norm (i.e. the norm of
	the gradient scaled by the size of the vector x and by
	f(x)) is less than this value. (1.0e-2) 

    Scalar StepReductionFactor: Reduction factor for backtracking 
        line search. (0.5) 

    Scalar FcnDecreaseTol: Ratio.of actual to predicted reduction
        at which step will be accepted (Goldstein-Armijo condition) (1.e-4) 

    Scalar SlopeDecreaseTol: Slope at approximate minimizer must be positive 
	or be reduced by a factor of SlopeDecreaseTol. (0.9)

    int MaxSample: Maximum number of function evaluations to perform in 
        line search (8)

    Scalar FractionOfMaxStep: Fraction of step to boundary of feasible set used
	to define maximum permissible step. (0.9) 

    Scalar FirstStep: Default initial step for line search.(1.0) 

    Scalar MinStep: Minimum allowable step.  If the algorithm takes a
	(relative) step less than this value in norm, it halts and reports
	"Possible convergence". (numeric_limits<Scalar>::min())

    Scalar MaxStep: Maximum allowable step.  If the algorithm takes too many
	consecutive steps of length MaxStep, it assumes that the
	iterates are diverging and reports "Possible divergence". 
	(numeric_limits<Scalar>::max())

    int CscMaxLimit:Maximum number of steps of length MaxStep allowed before
	the algorithm decides that the iterates are diverging. (5)

    int DispFlag: Display level.  This determines how much information should
	be displayed during the execution of the algorithm. Possible
	values are: 
	0 - No output; 
	1 - Function value and gradient norm after final iteration; 
	2 - Function value and gradient norm after every iteration.
        (0)

    int DispPrecision: Display precision---the number of digits sent to the
	screen. (6)

    Scalar InvHessianScale: Scale for the initial inverse Hessian approximation.  The
	initial inverse Hessian approximation is a multiple of the
	identity.  This scale is adjusted after each iteration using
	information gained from the step taken; however, its initial
	value must be provided. (1.0)

    int MaxUpdates: Maximum number of BFGS updates. See paper by Liu and Nocedal
	for details. (4)

    Scalar ZeroDivideTol: Tolerance for divide-by-zero - default is ProtDiv 
        overflow test (0.0)

    Scalar BracketIncrease: factor by which to increase step length if initial step 
        is too small (2.0)
    
 */

#include "table.hh"
#include "except.hh"

// <ME?> some conflict with predefined min and max requires this...?
#undef min
#undef max


namespace RVLUmin {

  using RVL::Table;
  using RVL::RVLException;

  template<class Scalar>
  class UMinTable: public Table {

  public:

    ~UMinTable() {}

    UMinTable(string fname = "") 
      : Table(fname) {
      
      Scalar ftmp;
      int itmp;
      string stmp;

      try {
	if (this->getValue("UpdateAlg",stmp)) {
	  stmp="LBFGS";
	  this->putValue("UpdateAlg",stmp);
	}
	if (this->getValue("UpdateAlg",stmp)) {
	  stmp="TRCG";
	  this->putValue("UpdateAlg",stmp);
	}
	if (this->getValue("LineSearchAlg",stmp)) {
	  stmp="BT";
	  this->putValue("LineSearchAlg",stmp);
	}
	if (this->getValue("MaxItn",itmp)) {
	  itmp = 10000;
	  this->putValue("MaxItn",itmp);
	}
	if (this->getValue("Typf",ftmp)) {
	  ftmp = 1.0;
	  this->putValue("Typf",ftmp);
	}
	if (this->getValue("TypxNorm",ftmp)) {
	  ftmp = 1.0;
	  this->putValue("TypxNorm",ftmp);
	}
	if (this->getValue("MaxUpdates",itmp)) {
	  itmp = 4;
	  this->putValue("MaxUpdates",itmp);
	}
	if (this->getValue("MaxSample",itmp)) {
	  itmp=8;
	  this->putValue("MaxSample",itmp);
	}
	if (this->getValue("BracketIncrease",ftmp)) {
	  ftmp=2.0;
	  this->putValue("BracketIncrease",ftmp);
	}
	if (this->getValue("GradTol",ftmp)) {
	  ftmp = 1.0e-2;
	  this->putValue("GradTol",ftmp);
	}
	if (this->getValue("FirstStep",ftmp)) {
	  ftmp = 1.0;
	  this->putValue("FirstStep",ftmp);
	}
	if (this->getValue("MinStep",ftmp)) {
	  this->putValue("MinStep", numeric_limits<Scalar>::min());
	}
	if (this->getValue("MaxStep",ftmp)) {
	  this->putValue("MaxStep", numeric_limits<Scalar>::max());
	}
	if (this->getValue("CscMaxLimit",itmp)) {
	  itmp = 5;
	  this->putValue("CscMaxLimit",itmp);
	}
	if (this->getValue("StepReductionFactor",ftmp)) {
	  ftmp=0.5;
	  this->putValue("StepReductionFactor",ftmp);
	}
	if (this->getValue("FcnDecreaseTol",ftmp)) {
	  ftmp = 1.e-4;
	  this->putValue("FcnDecreaseTol",ftmp);
	}
	if (this->getValue("SlopeDecreaseTol",ftmp)) {
	  ftmp = 0.9;
	  this->putValue("SlopeDecreaseTol",ftmp);
	}
	if (this->getValue("FractionOfMaxStep",ftmp)) {
	  ftmp = 0.9;
	  this->putValue("FractionOfMaxStep",ftmp);
	}
	if (this->getValue("DispFlag",itmp)) {
	  itmp = 3;
	  this->putValue("DispFlag",itmp);
	}
	if (this->getValue("DispPrecision",itmp)) {
	  itmp = 6;
	  this->putValue("DispPrecision",itmp);
	}
	if (this->getValue("InvHessianScale",ftmp)) {
	  ftmp = 1.0;
	  this->putValue("InvHessianScale",ftmp);
	}
	if (this->getValue("ZeroDivideTol",ftmp)) {
	  ftmp = 0.0;
	  this->putValue("ZeroDivideTol",ftmp);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from UMinTable constructor\n";
	throw e;
      }
    }

    UMinTable(const UMinTable<Scalar> & t)
      : Table(t) {}
   
    
  };

}

#endif
