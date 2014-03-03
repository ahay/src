/*************************************************************************

Copyright Rice University, 2011.
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

#ifndef __RVLALG_LBFGSBT_H
#define __RVLALG_LBFGSBT_H

#include "table.hh"
#include "lbfgsalg.hh"
#include "lnsrchBT.hh"
#include "uminstep.hh"


namespace RVLUmin {

  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::Functional;
  using RVL::FunctionalEvaluation;
  using RVL::Table;
  using RVL::RVLException;

/** Limited memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS)
    quasi-Newton optimization with geometric backtracking line search
    globalization

    For details of the LBFGS algorithm, see the paper
    
    "Updating Quasi-Newton Matrices with Limited Storage" by Jorge Nocedal,
    Math. of Computation, Vol. 35, no. 151, p.p. 773--782.

    Approximates solution to unconstrained optimization problem \f[
    \min_{x} f(x) \f] via generalized secant update \f[ x_+ = x_c -
    \alpha B \nabla f(x_c). \f]  \f$B\f$ is a limited-memory secant
    approximation to the inverse Hessian of \f$f\f$, parametrized by
    maximum rank. \f$\alpha\f$ is a step length.

    Note that this problem setting makes sense only for functions
    \f$f\f$ taking real values. Accordingly, the implementation
    includes a compile-time check that the Scalar template parameter
    designates a real field (i.e. double or float).

    This implementation globalizes convergence to a local min via a
    geometric backtracking linesearch method, impemented in
    RVLUmin::BacktrackingLineSearchAlg. Line search methods are
    described in Nocedal and Wright, <i>Numerical Optimization</i>,
    Springer 1999. For additional documentation, see
    RVLUmin::BacktrackingLineSearchAlg docs.

    Structure and function: Data members of three types combine to implement LBFGS:
<ul>
<li>RVLUmin::LBFGSDir - search direction update, combines update of inverse Hessian approximation (\f$B\f$, above) with computation of the BFGS search direction \f$B\nabla f(x_c)\f$</li>
<li>RVLUmin::BacktrackingLineSearchAlg - determines step length \f$\alpha\f$</li>
<li>RVLAlg::GradientThresholdIterationTable - monitors number of iterations and norm of gradient, stops iteration if max iteration count exceeded, or if gradient norm falls below threshold (absolute gradient norm tolerance), or if gradient norm falls below threshold relative to initial value at start of iteration (relative gradient norm tolerance)</li>
</ul>

    Usage: initialize solution vector (constructor argument x)
    externally; construct LBFGS object; call LBFGSBT::run(). On
    return, x has been updated to estimated solution.

    Parameters for LBFGS iteration:
    <ul>
    <li> f - function to be minimized (RVL::Functional)</li>
    <li> x - solution RVL::Vector - initial guess on call, estimated solution on return</li>
    <li> _ihs - inverse Hessian scale - overall scale factor, so initial Hessian is this Scalar multiple of identity operator. Typical value = 1.0</li>
    <li> _maxits - max permitted number of LBFGS iterations. Typical value = 20</li>
    <li> _mud - max stored BFGS updates - inverse Hessian approximation has this rank (at most). Typical value = 5.</li>
    <li> _agradtol - stopping tolerance for gradient norm, absolute. Typical value = 0.0</li>
    <li> _rgradtol - stopping tolerance for gradient norm, relative to initial gradient norm. Typical value = 0.01</li>
    <li> _str - verbose output stream
    </ul>

    NOTES: 
    <ul>

    <li> Setting _mud = 0 results in rank zero update, meaning that
    the Hessian is a multiple of the identity (the multiple being the
    inverse Hessian scale parameter, _ihs). That is, with this choice,
    the algorithm executes the Steepest Descent method, with
    backtracking line search globalization</li>

    <li> The inverse hessian scale parameter (_ihs) simply scales the
    line search. We include it to conform to common BFGS
    implementations, however a reasonable choice appears to be
    1.0.</li>

    <li> Setting either stopping tolerance to 0.0 disables it. To stop
    entirely based on the gradient norm reduction relative to the
    initial iterate, set _agradtol to 0.0.</li>
    </ul>

    Parameters for backtracking line search (for detailed description
    and notes, see docs for RVLUmin::BacktrackingLineSearchAlg):

    <ul>
    <li> _maxsamp - max number of steps permitted in each line search

    <li> _disp - verbosity flag - false = no output, true = function
    value, gradient norm at each iteration, report of line search

    <li> _sl1 - first line search step
    <li> _eta1 - lower G-A parameter
    <li> _eta2 - upper G-A parameter
    <li> _gamma1 - line search backtrack factor
    <li> _gamma2 - line search extrapolation factor ("internal doubling")
    <li> _maxfrac - fraction of max step to boundary permitted
    </ul>

    Typical use case: see <a href="../../testsrc/testlbfgs.cc">functional test source</a>.
 */

  template<typename Scalar>
  class LBFGSBT: public Algorithm {
    
  private:

    // parameters for LBFGSDir
    Scalar ihs;        // inverse Hessian scale
    int mud;           // max updates

    // parameters for BT line search
    int maxsamp;       // max function evaluations
    bool disp;         // display flag 
    Scalar sl1;        // length of first step
    Scalar minsteptol; // minimum permitted step length (fraction of prev step)
    Scalar eta1;       // First GA scale: min acceptable decrease
    Scalar eta2;       // Second GA scale: good decrease
    Scalar gamma1;     // First BT factor: shrink step if decrease not acceptable
    Scalar gamma2;     // Second BT factor: expand step if decrease good
    Scalar maxfrac;    // fraction of max step to attempt

    // parameters for loop alg
    int maxits;        // max number of LBFGS steps
    Scalar agradtol;   // absolute gradient stopping threshhold
    Scalar rgradtol;   // relative gradient stopping threshhold (to initial)
    // also uses disp
    ostream & str;
    FunctionalEvaluation<Scalar> fx;
    LBFGSDir<Scalar> dir;
    BacktrackingLineSearchAlg<Scalar> ls;
    GradientThresholdIterationTable<Scalar> ctbl;

    LBFGSBT();
    LBFGSBT(const LBFGSBT<Scalar> &);

  public:

    /** constructor, first form:

	parameters:
	@param f - function to be minimized (RVL::Functional)
	@param x - solution RVL::Vector - initial guess on call, estimated solution on return
	@param _ihs - inverse Hessian scale - overall scale factor, so initial Hessian is this Scalar multiple of identity operator
	@param _maxits - max number of LBFGS iterations
	@param _mud - max stored BFGS updates - stored inverse Hessian approximation has this rank (at most)
	@param _maxsamp - max number of steps permitted in each line search
	@param _disp - verbosity flag - false = no output, true = function value, gradient norm at each iteration, report of line search
	@param _sl1 - first line search step
	@param _eta1 - lower G-A parameter
	@param _eta2 - upper G-A parameter
	@param _gamma1 - line search backtrack factor
	@param _gamma2 - line search extrapolation factor ("internal doubling")
	@param _maxfrac - fraction of max step to boundary permitted
	@param _agradtol - stopping tolerance for gradient norm, absolute
	@param _rgradtol - stopping tolerance for gradient norm, relative to initial gradient norm
	@param _str - verbose output unit

    */

    LBFGSBT(Functional<Scalar> const & f,
	    Vector<Scalar> & x,
	    Scalar _ihs,
	    int _mud,
	    int _maxsamp,
	    bool _disp,
	    Scalar _sl1,
	    Scalar _eta1,
	    Scalar _eta2,
	    Scalar _gamma1,
	    Scalar _gamma2,
	    Scalar _maxfrac,
	    Scalar _minsteptol,
	    int _maxits,
	    Scalar _agradtol,
	    Scalar _rgradtol,
	    ostream & _str = cout) 
      : ihs(_ihs), mud(_mud),
	maxsamp(_maxsamp),
	disp(_disp),
	sl1(_sl1),
	eta1(_eta1),
	eta2(_eta2),
	gamma1(_gamma1),
	gamma2(_gamma2),
	maxfrac(_maxfrac),
	minsteptol(_minsteptol),
	maxits(_maxits),
	agradtol(_agradtol),
	rgradtol(_rgradtol),
	str(_str),
	fx(f,x),
	dir(fx.getDomain(),ihs,mud,str),
	ls(maxsamp,disp,sl1,minsteptol,eta1,eta2,gamma1,gamma2,maxfrac,str),
	ctbl(fx,maxits,agradtol,rgradtol,str)
    { testRealOnly<Scalar>(); }

    /** constructor, second form - scalar arguments passed in Table. */
    LBFGSBT(Functional<Scalar> const & f,
	    Vector<Scalar> & x,
	    Table const & t,
	    ostream & _str = cout)
      : ihs(getValueFromTable<Scalar>(t,"BFGS_InvHessianScale")),
	mud(getValueFromTable<int>(t,"BFGS_MaxUpdates")),
	maxsamp(getValueFromTable<int>(t,"LS_MaxSample")),
	disp(getValueFromTable<bool>(t,"DispFlag")),
	sl1(getValueFromTable<Scalar>(t,"LS_FirstStep")),
	minsteptol(getValueFromTable<Scalar>(t,"LS_MinStepTol")),
	eta1(getValueFromTable<Scalar>(t,"MinDecrease")),
	eta2(getValueFromTable<Scalar>(t,"GoodDecrease")),
	gamma1(getValueFromTable<Scalar>(t,"StepDecrFactor")),
	gamma2(getValueFromTable<Scalar>(t,"StepIncrFactor")),
	maxfrac(getValueFromTable<Scalar>(t,"LS_FractionOfMaxStep")),      
	maxits(getValueFromTable<int>(t,"MaxItn")),
	agradtol(getValueFromTable<Scalar>(t,"AbsGradTol")),
	rgradtol(getValueFromTable<Scalar>(t,"RelGradTol")),
	str(_str),
	fx(f,x),
	dir(fx.getDomain(),ihs,mud,str),
	ls(maxsamp,disp,sl1,minsteptol,eta1,eta2,gamma1,gamma2,maxfrac,str),
	ctbl(fx,maxits,agradtol,rgradtol,str)
    { testRealOnly<Scalar>(); }
    int getCount() { return ctbl.getCount(); }

    void run() {
      try {
	// put together alg, terminator, and loopalg.
	// note that BFGS update rule is also a terminator
	UMinStepLS<Scalar> step(fx,dir,ls,str);
	OrTerminator stop(step,ctbl);
	LoopAlg Umin(step, stop);
	Umin.run();
      }

      catch (RVLException & e) {
	e<<"\ncalled from LBFGSBT constructor\n";
	throw e;
      }
    }

    ~LBFGSBT() {}

    /** supplied to provide access to any intermediate data that a subclass
	of Functional may make available. Can extract const reference to current
	Functional, as constructed by FunctionalEvaluation, via the getFcnl() method.
	A cast will be required to extract any further subclass attributes.
    */
    FunctionalEvaluation<Scalar> const & getFunctionalEvaluation() const { return fx; }
    
  };

}
#endif
