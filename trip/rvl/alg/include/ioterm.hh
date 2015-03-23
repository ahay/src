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

#ifndef __RVLALG_IO_TERMINATOR
#define __RVLALG_IO_TERMINATOR

#include "std_cpp_includes.hh"
#include "alg.hh"
#include "scalarterm.hh"
#include "functional.hh"

/**  ioterm.H
     All Terminators which perform some I/O when queried.
     This includes terminators which prompt the user for some
     action as well as those which print out messages and tables.
*/

namespace RVLAlg {
  using namespace RVL;
/** Prints a message to the output stream and reads a single character response 
    from the input stream

    returns true if  reponse == 'y' or 'Y'
            false    otherwise
*/


class IOTerminator: public Terminator {
public:
  
  IOTerminator() : ins(cin), outs(cout) { 
    s = new char[30];  
    strncpy(s, "Do you wish to stop?", 30); 
  }

  IOTerminator(char t[]) : ins(cin), outs(cout) { 
    int l = strlen(t);
    s = new char[l];  
    strncpy(s, t, l);
  }

  IOTerminator(istream & in_, ostream & out_) : ins(in_), outs(out_) 
  {  
    s = new char[30];  
    strncpy(s, "Do you wish to stop iterating?", 30); 
  }

  IOTerminator(char t[], istream & in_, ostream & out_) : ins(in_), outs(out_) 
  { 
    int l = strlen(t);
    s = new char[l];  
    strncpy(s, t, l);
  }

  virtual ~IOTerminator() {
    if( s != NULL )
      delete [] s;
  }
    
  virtual bool query() {
    char response;
    outs << s << " y/n : ";
    ins >> response;
    ins.ignore(10, '\n');     // Clear buffer? 

    if ((response == 'y')||(response == 'Y'))
      return 1;
    else
      return 0;
  }

protected:
  char * s;
  istream & ins;
  ostream & outs;

};

/** Rather odd terminator which does nothing except 
    output the value of a vector when queried
    
    NOTE: always returns 0;
*/

template< class Scalar>
class VecWatchTerminator: public Terminator {
public:
  VecWatchTerminator( Vector<Scalar> & x ) 
    : x_(x), outs_(cout){}

  VecWatchTerminator( Vector<Scalar> & x, ostream & outs ) 
    : x_(x), outs_(outs){}

  virtual bool query() {
    x_.write(outs_);
    return 0;
  }

private: 
  Vector<Scalar> & x_;
  ostream & outs_;
};

  /** This terminator never says to stop, but simply
      prints the current functional value.
  */
template< class Scalar >
class IterationTable: public Terminator {
public:
  IterationTable( FunctionalEvaluation<Scalar> & _fx, 
		  ostream & _out = cout) 
    : count(0), fx(_fx), out(_out)
  {
    out << setw(10) << "Iteration " << "|" 
	 << setw(13) <<"f(x)" << endl;
  }

  bool query() {
    out << setw(10) << count << ' '; 
    out << setw(13) << setprecision(8) << fx.getValue() << endl;
    count++;
    return false;
  }
    
protected:
  int count;
  FunctionalEvaluation<Scalar> & fx;
  ostream & out;

};

  /** This terminator never says to stop, but simply
      prints the current functional value and 
      an associated scalar (ussually the step size)
  */
template< class Scalar >
class SteppedIterationTable: public Terminator {
public:
  SteppedIterationTable(FunctionalEvaluation<Scalar> & _fx,
			Scalar & _step,
			ostream & _out = cout) 
    : count(0), fx(_fx), step(_step), out(_out)
  {
    out << setw(10) << "Iteration " << "|" 
	 << setw(13) <<"f(x)" << "|"
	 << setw(13) << "Delta" << endl;
  }

  bool query() {
    
    Scalar fxval;

    fxval = fx.getValue();

    out << setw(10) << count << ' '; 
    out << setw(13) << setprecision(8) << fxval << ' ';
    out << setw(13) << setprecision(8) << step << endl;
    count++;
    return false;
  }
    
protected:
  int count;
  FunctionalEvaluation<Scalar> & fx;
  Scalar & step;
  ostream & out;

};

  /** This terminator behaves like a combined CountTerminator and
      NormGradientTerminator, with the added side effect of printing
      an iteration table with the current functional value and norm of
      the gradient.  Tests both absolute and relative (to initial
      value) gradient norm tolerances, returns true when either
      descends below given tolerance. To force absolute gradient norm
      to descend below atol, set rtol=1; to force relative gradient
      norm to descend below rtol, set atol=0. Defaults are set to be
      satisfied at first iteration, for safety's sake.

      Template on Scalar type, with absolute value type = abstype

      @param[in] fx - (FunctionalEvaluation)
      @param[in] maxcount - (int) max number of iterations permitted; default = 1
      @param[in] atol - (abstype) absolute gradient tolerance; default = max value of abstype
      @param[in] rtol - (abstype) relative gradient tolerance; default = abstype one
      @param[in] doit - (bool) print table headings; default = true
      @param[in] out - (ostream) output stream; default = cout
       
  */
template< class Scalar >
class GradientThresholdIterationTable: public CountTerminator {
  
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
  
public:
  GradientThresholdIterationTable(FunctionalEvaluation<Scalar> & _fx, 
				  int maxcount = 1, 
				  atype _atol = numeric_limits<atype>::max(),
				  atype _rtol = ScalarFieldTraits<atype>::One(),
				  ostream & _out = cout,
				  bool _doit = true) 
      : CountTerminator(maxcount), atol(_atol), rtol(_rtol), init(true), ngfx0(ScalarFieldTraits<atype>::Zero()), ingfx(ScalarFieldTraits<atype>::One()), doit(_doit),  fx(_fx), out(_out) {
    if (doit) {
      out << setw(10) << "Iteration " << "|" 
	  << setw(13) <<"f(x)   " << "|" 
	  << setw(13) << "||gradf(x)||" << endl;
    }
  }

  bool query() {
    
    Scalar fxval;
    atype ngfx;
    bool stop = CountTerminator::query();
    
    fxval = fx.getValue();
    ngfx  = fx.getGradient().norm();

    if (init) {
      if (ProtectedDivision<atype>(ScalarFieldTraits<atype>::One(),ngfx,ingfx)) {
	out<<"Initial gradient below roundoff - success!! (?)\n";
	return true;
      }
      ngfx0=ngfx;
      fxval0=fxval;
      init=false;
    }
    
    bool astop = (ngfx <= atol);
    bool rstop = (ngfx*ingfx <= rtol );
    stop = (stop || astop || rstop );
    
    if (doit) {
      if( !stop ) { //print out current results
	out << scientific << setw(10) << getCount() << ' '; 
	out << scientific << setw(13) << setprecision(8) << fxval << ' '
	    << scientific << setw(13) << setprecision(8) << ngfx << endl;
      }
      else { // print out hline, then final results 
	out << "-----------------------------------------------" << endl;
	if (astop || rstop) {
	  out << "Convergence: gradient norm below tolerance after "<< setw(10) << getCount() << " iterations" <<endl;
	}
	else {
	  out <<"Termination: maximum number of iterations exceeded"<<endl;
	}
	out << "Function Value           = " << scientific << setw(13) << setprecision(8) << fxval <<endl;
	out << "Initial Function Value   = " << scientific << setw(13) << setprecision(8) << fxval0 <<endl;
	atype frat;
	if ((fxval >= ScalarFieldTraits<atype>::Zero()) && 
	    (fxval0 > ScalarFieldTraits<atype>::Zero()) &&
	    (!(ProtectedDivision<atype>(fxval,fxval0,frat)))) {
	  out << "Function Value Reduction = " << scientific << setw(13) << setprecision(8) << frat <<endl;
	}
	out << "Gradient Norm            = " << scientific << setw(13) << setprecision(8) << ngfx << endl;
	out << "Initial Gradient Norm    = " << scientific << setw(13) << setprecision(8) << ngfx0 << endl;
	atype grat;
	if ((ngfx >= ScalarFieldTraits<atype>::Zero()) && 
	    (ngfx0 > ScalarFieldTraits<atype>::Zero()) &&
	    (!(ProtectedDivision<atype>(ngfx,ngfx0,grat)))) {
	  out << "Gradient Norm Reduction  = " << scientific << setw(13) << setprecision(8) << grat <<endl;
	}
	out << "Absolute Tolerance       = " << scientific << setw(13) << setprecision(8) << atol << endl;
	out << "Relative Tolerance       = " << scientific << setw(13) << setprecision(8) << rtol << endl;
      }
    }

    return stop;
  } 
      
protected:
  atype atol;
  atype rtol;
  mutable bool init;
  mutable atype fxval0;
  mutable atype ngfx0;
  mutable atype ingfx;
  
  bool doit;
  FunctionalEvaluation<Scalar> & fx;
  ostream & out;
  
};

 /** This terminator behaves like a combined CountTerminator
      and MinTerminator, with the added side
      effect of printing an iteration table with the
      value of the watched Scalar.

      Returns stop when the count exceeds maxcount or the Scalar
      falls below the specified tolerance.
  */
template< class Scalar >
class CountingThresholdIterationTable: public CountTerminator {
public:

  typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  CountingThresholdIterationTable(int maxcount,
                                  const atype & val,
                                  atype tol,
                                  string & scalarname,
                                  ostream & out = cout) 
    : tol_(tol), val_(val), out_(out), CountTerminator(maxcount)  
  {
    out_ << setw(10) << "Iteration " << " | " 
         << setw(scalarname.size()+5) << scalarname << endl;
  }

  CountingThresholdIterationTable(int maxcount,
                                  const atype & val,
                                  atype tol,
                                  ostream & out = cout) 
    : tol_(tol), val_(val), out_(out), CountTerminator(maxcount)  
  {
    out_ << setw(10) << "Iteration " << "|" 
         << setw(17) << "Scalar Value" << endl;
  }

  bool query() {
    
    bool stop = CountTerminator::query();
    if( !stop ) { //print out current results
      out_ << setw(10) << getCount() << ' '; 
      out_ << setw(17) << setprecision(8) << val_ << endl;
    }
    else { // print out hline, then final results 
      out_ << "-----------------------------------------------" << endl;
      out_ << setw(10) << getCount() << ' ' 
           << setw(17) << setprecision(8) << val_ << endl;
    }
    return (stop || (val_ <= tol_ ));
  }
    
    
protected:
  atype tol_;
  const atype & val_;
  ostream & out_;

};



 /** Vector version of CountingThreshholdIterationTable Note separate
     initialization - allows object to be instantiated in member
     initialization list, while names, numbers, etc. are assigned
     later, eg. in constructor body
  */
template< class Scalar >
class VectorCountingThresholdIterationTable: public CountTerminator {
public:

  typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  VectorCountingThresholdIterationTable(int maxcount,
					vector<string> const & _names,
					vector<atype *> const & _nums,
					vector<atype> const & _tols,
					ostream & _out = cout) 
    :  CountTerminator(maxcount), tols(_tols), nums(_nums), names(_names), out(_out)
  {
    if ((nums.size() != names.size()) ||
	(tols.size() != names.size())) {
      RVLException e;
      e<<"Error: VectorCountingThresholdIterationTable constructor\n";
      e<<"input name, value, and tol arrays must have same length\n";
      throw e;
    }
  }

  void init() {
    out << setw(10) << "Iteration  ";
    //    cerr << setw(10) << "Iteration  ";
    for (size_t i=0;i<names.size();i++) {
      out  << " |  "<< setw(names[i].size()) << names[i];
      //      cerr  << " |  "<< setw(names[i].size()+1) << names[i];
    }
    out << endl;
    //    cerr << endl;
  }

  bool query() {
    
    bool stop = CountTerminator::query();

    if( !stop ) { //print out current results
      out << setw(10) << getCount() << ' ';
      for (size_t i=0;i<names.size();i++)  
	out << setw(4+names[i].size()) << setprecision(8) << *(nums[i]);
      out << endl;
    }
    else { // print out hline, then final results 
      out << "-----------------------------------------------" << endl;
      out << setw(10) << getCount() << ' ' ;
      for (size_t i=0;i<names.size();i++)  
	out << setw(4+names[i].size()) << setprecision(8) << *(nums[i]);
      out << endl;
    }
    for (size_t i=0; i<names.size();i++) {
      if (*(nums[i]) < tols[i]) out<<names[i]<<" = "<<*(nums[i])<<" below tol = "<<tols[i]<<endl;
      stop = stop || (*(nums[i]) < tols[i] );
    }
    return stop;
  }
    
    
protected:
  vector<atype> tols;
  vector<atype *> nums;
  vector<string> names;
  ostream & out;

private:

  VectorCountingThresholdIterationTable();
  VectorCountingThresholdIterationTable( VectorCountingThresholdIterationTable const & );
};

/** This terminator behaves like a combined CountTerminator
      and NormTerminator, with the added side
      effect of printing an iteration table with the
      current value of the norm of the vector.
  */
template< class Scalar >
class CountingNormTable: public CountTerminator {
public:
  typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  CountingNormTable(int maxcount, 
		    const Vector<Scalar> & _x, 
		    NormRetType _tol,
		    bool _doit = true,
		    ostream & _out = cout) 
    : x(_x), tol(_tol), doit(_doit), out(_out), CountTerminator(maxcount) {
    if (doit) {
      out << setw(10) << "Iteration " << "|" 
	  << setw(13) << "||x||" << endl;
    }
  }

  bool query() {
    
    NormRetType nx;
    bool stop = CountTerminator::query();

    nx = x.norm();

    if (doit) {
      if( !stop ) { //print out current results
	out << setw(10) << getCount() << ' '; 
	out << setw(13) << setprecision(8) << nx << endl;
      }
      else { // print out hline, then final results 
	out << "-----------------------------------------------" << endl;
	out << setw(10) << getCount() << ' '; 
	out << setw(13) << setprecision(8) << nx << endl;
      }
    }
    return (stop || (nx <= tol ));
  }
    
    
protected:
  bool doit;
  ostream & out;
  const Vector<Scalar> & x;
  NormRetType tol;
};

}

#endif
