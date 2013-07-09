/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007, 2008
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
#ifndef __TSOPT_OP
#define __TSOPT_OP

#include "jet.hh"
#include "op.hh"

namespace TSOpt {

  using RVL::Space;
  using RVL::Vector;
  using RVL::FunctionObject;
  using RVL::FunctionObjectConstEval;
  using RVL::Operator;
  using RVL::Writeable;


  /** Adaptor class, making an RVL Operator out of TSOpt's jet class.
      Uses jet object to define an "apply" operator, "applyDeriv" 
      operator and "applyAdjDeriv" operator. 
  */
  template
  <
    typename Scalar, 
    typename State,
    class InitCreationPolicy,
    class FinalCreationPolicy
  >

  class SimOp: public Operator<Scalar>, 
	       public InitCreationPolicy, 
	       public FinalCreationPolicy 
  {
      
  private:
    
    Sim1Jet<State> & j;
    Space<Scalar> const & dom;
    Space<Scalar> const & rng;
    
    SimOp();
      
  protected:
      
    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {

      //      cout << "In SimOp::apply" << endl;
      try {
	if (this->getDomain() != x.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::apply\n";
	  e<<"input vector not in domain\n";
	  e<<"input vector:\n";
	  x.write(e);
	  e<<"domain\n";
	  this->getDomain().write(e);
	}
	if (this->getRange() != y.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::apply\n";
	  e<<"output vector not in range\n";
	  e<<"output vector:\n";
	  y.write(e);
	  e<<"range\n";
	  this->getRange().write(e);
	}
	
	FunctionObjectConstEval * f = InitCreationPolicy::create(j.getFwd().getState());
	x.eval(*f);
	delete f;
	j.getFwd().run();
	FunctionObject * g = FinalCreationPolicy::create(j.getFwd().getState());
	y.eval(*g);
	delete g;
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::SimOp::apply\n";
	throw e;
      }

      //      cout << "in SimOp::apply" << endl;
      //      cout << "x = " << endl;
      //      x.write(cout)

    }
      
    void applyDeriv(const Vector<Scalar> & x, 
		    const Vector<Scalar> & dx,
		    Vector<Scalar> & dy) const {

      //      cout << "in SimOp::applyDer" << endl;

      try {
	if (this->getDomain() != x.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyDeriv\n";
	  e<<"input vector not in domain\n";
	  e<<"input vector:\n";
	  x.write(e);
	  e<<"domain\n";
	  this->getDomain().write(e);
	}

	if (this->getDomain() != dx.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyDeriv\n";
	  e<<"input vector not in domain\n";
	  e<<"input vector:\n";
	  dx.write(e);
	  e<<"domain\n";
	  this->getDomain().write(e);
	}

	if (this->getRange() != dy.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyDeriv\n";
	  e<<"output vector not in range\n";
	  e<<"output vector:\n";
	  dy.write(e);
	  e<<"range\n";
	  this->getRange().write(e);
	}


	
	FunctionObjectConstEval * f = InitCreationPolicy::create(j.getFwd().getState());
	x.eval(*f);
	delete f;

	//	cout << "B: linState = " << endl;
	//	j.getLin().getState().write(cout);

	
	FunctionObjectConstEval * df = InitCreationPolicy::create(j.getLin().getState());
	dx.eval(*df);

	//	cout << "A: linState = " << endl;
	//	j.getLin().getState().write(cout);


	//	cout << "B: dx = " << endl;
	//	dx.write(cout);

	delete df;

	j.getLin().run();
	
	//	cout << "A: dx = " << endl;
	//	dx.write(cout);
      
	FunctionObject * g = FinalCreationPolicy::create(j.getLin().getState());
	dy.eval(*g);
	delete g;

      }

      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::SimOp::applyDeriv\n";
	throw e;
      }

      /*
    
      cout << "x = " << endl;
      x.write(cout);
      cout << "dx = " << endl;
      dx.write(cout);
      cout << "dy = " << endl;
      dy.write(cout);
      */
    }
      
    void applyAdjDeriv(const Vector<Scalar> & x, 
		       const Vector<Scalar> & dy,
		       Vector<Scalar> & dx) const {

      //      cout << "in SimOp::applyAdjDeriv" << endl;
      if (this->getDomain() != x.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyAdjDeriv\n";
	  e<<"input vector not in domain\n";
	  e<<"input vector:\n";
	  x.write(e);
	  e<<"domain\n";
	  this->getDomain().write(e);
	}
 
	if (this->getDomain() != dy.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyAdjDeriv\n";
	  e<<"input vector not in domain\n";
	  e<<"input vector:\n";
	  dy.write(e);
	  e<<"domain\n";
	  this->getDomain().write(e);
	}

	if (this->getRange() != dx.getSpace()) {
	  RVLException e;
	  e<<"Error: TSOpt::SimOp::applyAdjDeriv\n";
	  e<<"output vector not in range\n";

	  e<<"output vector:\n";
	  dy.write(e);
	  e<<"range\n";
	  this->getRange().write(e);
	}

	FunctionObjectConstEval * f = InitCreationPolicy::create(j.getFwd().getState());
	x.eval(*f);
	delete f;

	FunctionObjectConstEval * df = InitCreationPolicy::create(j.getAdj().getState());
	dy.eval(*df);
	delete df;

	j.getFwd().run();
      

	j.getAdj().setTime(j.getFwd().getTime());
       

	j.getAdj().run();

	FunctionObject * g = FinalCreationPolicy::create(j.getAdj().getState());
	dx.eval(*g);
	delete g;

	/*

	cout << "x = " << endl;
	x.write(cout);
	cout << "dy = " << endl;
	dy.write(cout);
	cout << "dx = " << endl;
	dx.write(cout);
	*/


    }
      
    Operator<Scalar> * clone() const { 
      return new SimOp<Scalar,State,InitCreationPolicy,FinalCreationPolicy>(*this);
    }
      
  public:
      
    SimOp(Space<Scalar> const & _dom,
	  Space<Scalar> const & _rng,
	  Sim1Jet<State> & _j)
      : dom(_dom), rng(_rng), j(_j) {}
    SimOp(SimOp
	  <
	  Scalar,
	  State,
	  InitCreationPolicy,
	  FinalCreationPolicy 
	  > const & x)
      : dom(x.dom), rng(x.rng), j(x.j) {}
    virtual ~SimOp() {}
      
    virtual const Space<Scalar> & getDomain() const { return dom; }
    virtual const Space<Scalar> & getRange() const { return rng; }
      
    ostream & write(ostream & str) const {
      return str;
    }
  };
}


#endif
