/****************************************
 added by Yin Huang at May 25, 2013 
 Test LinFitLS with nonlinear operator cube
*************************************************************************

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


#include "rnmat.hh"
#include "adjtest.hh"
#include "cgnealg.hh"
#include "LinFitLS.hh"
#include "chebalg.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;


template<class Scalar>
class cube: public Operator<Scalar>{
private:
    RnSpace<Scalar> Dom;
    RnSpace<Scalar> Rng;
    mutable bool applied;
    cube(){}
public:
    cube(RnSpace<Scalar> _dom, RnSpace<Scalar> _rng,bool _applied=false): Dom(_dom), Rng(_rng), applied(_applied){}
    cube(cube const & x):applied(x.applied){};
    ~cube(){};
    
protected:
    void apply(const Vector<Scalar> & x,
               Vector<Scalar> & y) const{
        try {
            LocalVector<Scalar> vin(x);
            LocalVector<Scalar> vout(y);

            for (int i=0;i<vin.getSize();i++)
                vout.getData()[i] = (vin.getData()[i])*(vin.getData()[i])*(vin.getData()[i]);
            
            y.copy(vout);
        } catch (RVLException & e) {
            e<<"\ncalled from cube::apply\n";
            throw e;
        }
        
    };
    
    void applyDeriv(const Vector<Scalar> & x,
                    const Vector<Scalar> & dx,
                    Vector<Scalar> & dy) const{
        try {
            LocalVector<Scalar> in(x);
            LocalVector<Scalar> din(dx);
            LocalVector<Scalar> dout(dy);
            for (int i=0;i<in.getSize();i++)
                dout.getData()[i] = 3.0*(in.getData()[i])*(in.getData()[i])*(din.getData()[i]);
            dy.copy(dout);
        } catch (RVLException & e) {
            e<<"\ncalled from cube::applyDeriv\n";
            throw e;
        }
    };
    
    void applyAdjDeriv(const Vector<Scalar> & x,
                       const Vector<Scalar> & dy,
                       Vector<Scalar> & dx) const{
        try {
            LocalVector<Scalar> in(x);
            LocalVector<Scalar> dout(dx);
            LocalVector<Scalar> din(dy);
            for (int i=0;i<in.getSize();i++)
                dout.getData()[i] = 3.0*(in.getData()[i])*(in.getData()[i])*(din.getData()[i]);
            dx.copy(dout);
        } catch (RVLException & e) {
            e<<"\ncalled from cube::applyAdjDeriv\n";
            throw e;
        }
        
    };
    
    void applyDeriv2(const Vector<Scalar> & x,
                     const Vector<Scalar> & dx0,
                     const Vector<Scalar> & dx1,
                     Vector<Scalar> & dy) const{
        try {
            LocalVector<Scalar> in(x);
            LocalVector<Scalar> din0(dx0);
            LocalVector<Scalar> din1(dx1);
            LocalVector<Scalar> dout(dy);
            for (int i=0;i<in.getSize();i++)
                dout.getData()[i] = 6.0*(in.getData()[i])*(din0.getData()[i])*(din1.getData()[i]);
            dy.copy(dout);

        } catch (RVLException & e) {
            e<<"\ncalled from cube::applyDeriv2\n";
            throw e;
        }
        
    }
    
    void applyAdjDeriv2(const Vector<Scalar> & x,
                        const Vector<Scalar> & dx0,
                        const Vector<Scalar> & dy,
                        Vector<Scalar> & dx1) const{
        try {
            LocalVector<Scalar> in(x);
            LocalVector<Scalar> din0(dx0);
            LocalVector<Scalar> din1(dy);
            LocalVector<Scalar> dout(dx1);
            for (int i=0;i<in.getSize();i++)
                dout.getData()[i] = 6.0*(in.getData()[i])*(din0.getData()[i])*(din1.getData()[i]);
            dx1.copy(dout);
        } catch (RVLException & e) {
            e<<"\ncalled from cube::applyAdjDeriv2\n";
            throw e;
        }
        
    }
    const Space<Scalar> & getRange() const{return Rng;}
    
    const Space<Scalar> & getDomain() const{return Dom;}
    
    Operator<Scalar> * clone() const {return new cube<Scalar> (*this);}
    
    ostream & write(ostream & str) const{
        str<<"cube object\n";
        return str;
    }
    
};

template<typename T>
void LinFitTestNL(int rows, int cols, typename ScalarFieldTraits<T>::AbsType _nlvl,
		int maxit, 
		typename ScalarFieldTraits<T>::AbsType rtol, 
		typename ScalarFieldTraits<T>::AbsType nrtol, 
		typename ScalarFieldTraits<T>::AbsType Delta, 
		string testno, string datatype, string output, 
		ostream & str) {
  
  typedef typename ScalarFieldTraits<T>::AbsType atype;
  
  try {
    str<<endl;
    str<<"/****************************************************"<<endl;
    str<<" *            BEGIN UMin UNIT TEST "<<testno<<endl;
    str<<" * Purpose: test LinFitLS comp type ="<<datatype<<endl;
    str<<" *   noise level = "<<_nlvl<<endl;
    str<<" *   results reported in "<<output<<endl;
    str<<" ****************************************************/"<<endl; 
    str<<endl;

    // build output stream
    ofstream strout(output.c_str());
    strout<<scientific;

    str<<"1. Construct domain, range spaces: dom=rng"<<endl;
    RnSpace<T> dom(cols);
    RnSpace<T> rng(cols);
    str<<endl;

    str<<"2. Build operator cube (cube), do adjoint test\n";
    cube<T> nlop(dom,dom);
    RVLRandomize<T> rnd;

    str<<"3. Build sol, rhs vectors\n";
    Vector<T> sol(dom);
    Vector<T> x0(dom);
    Vector<T> x1(dom);
    Vector<T> est(dom);
    // initialize exact solution to random vector
    // x1 = x0 + dltx
    x0.eval(rnd);
    x1.eval(rnd);
    x1.scale(0.1);
    x1.linComb(1.0,x0);
      
    // Operator cube evaluate at x0 and x1
    OperatorEvaluation<T> Fx0(nlop,x0);
    OperatorEvaluation<T> Fx1(nlop,x1);
 
    // adj test for Deriv of Operator
    LinearOp<T> const & lop=Fx0.getDeriv();
      SymmetricBilinearOp<T> const & sblop = Fx0.getDeriv2();
      if (AdjointTest<T>(lop,rnd,strout)) {
          str<<"*** LinearOp adjoint test passed - result in "<<output<<"; proceed\n";
      }
      else {
          str<<"*** LinearOp adjoint test failed - result in "<<output<<"; exit\n";
          return;
      }
      str<<endl;

    // initialize estimated solution to zero
    est.zero();
    // create rhs, initialize to noisy version of exact sol image
    // rhs = F(x1)-F(x0)
    Vector<T> const & val0=Fx0.getValue();
      Vector<T> const & val1=Fx1.getValue();
      Vector<T> rhs(val1);
      rhs.linComb(-1.0,val0);
      
      Vector<T> err(dom);
      err.eval(rnd);
      T nlvl = _nlvl*rhs.norm()/err.norm();
      rhs.linComb(nlvl,err);
      str<<endl;
      
    str<<"4. Build, run CGNE algorithm \n";
    str<<endl;    
    strout<<endl<<"*******************************************************"<<endl;
    strout<<"* Test of LinFitLS Function"<<endl;
    strout<<"* max iterations       = "<<maxit<<endl;
    strout<<"* residual tolerance   = "<<rtol<<endl;
    strout<<"* normal res tolerance = "<<nrtol<<endl;
    strout<<"* trust radius         = "<<Delta<<endl;


    CGNEPolicyData<T> params;
    params.rtol=rtol;
    params.nrtol=nrtol;
    params.maxcount=maxit;
    params.verbose=false;

    // allocate rnorm, nrnorm scalars
    atype rnorm = numeric_limits<atype>::max();
    atype nrnorm = numeric_limits<atype>::max();

    // construct CGNEAlg
    CGNEPolicy<T> p;
    p.assign(params);
    CGNEAlg<T> * iter = p.build(est,lop,rhs,rnorm,nrnorm,strout);
    
    // run it
    iter->run();
    // compute gradient
    Vector<T> dltd(dom);
    Vector<T> grad(dom);
    lop.applyOp(est,dltd);
    dltd.linComb(-1.0,rhs);
    sblop.applyAdjOp(est,dltd,grad);
      

    // Generate identity preconditioner
    GenMat<T> B(dom,dom);
    for (int i=0;i<cols;i++)
          B.getElement(i,i)=ScalarFieldTraits<T>::One();
    // construct LinFitLS functional
    LinFitLS<T, CGNEPolicy<T>, CGNEPolicyData<T> > f(nlop,B,rhs,params,strout);

    FunctionalEvaluation<T> feval(f,x0);
    strout<<"\n*******************************************************";
    strout<<"\n Results from running CGNE\n";
    strout<<"value from LinFitLS    = "<<feval.getValue()<<endl;
    strout<<"rnorm from CGNE        = "<<rnorm<<endl;
    strout<<"rnorm from LinFitLS    = "<<sqrt(2.0*abs(feval.getValue()))<<endl;
    strout<<"difference rnorms      = "<<rnorm-sqrt(2.0*abs(feval.getValue()))<<endl;
      // compute correction term of grad
      CGNEPolicy<T> p0;
      p0.assign(params);
      est.zero();
      CGNEAlg<T> * iter0 = p0.build(est,lop,dltd,rnorm,nrnorm,strout);
      iter0->run();
      Vector<T> tmp0(dom);
      sblop.applyAdjOp(est,rhs,tmp0);
      grad.linComb(1.0, tmp0);
    strout << "Output Gradient Norm:\n";
    strout << "Grad norm from LinFitLS         = "<<feval.getGradient().norm() << endl;
    strout << "Grad norm computed in this test = "<<grad.norm() << endl;
    strout<<"*******************************************************"<<endl;

      
    str<< "5. build, run chebalg\n";
      ChebPolicyData<T> param1;
      param1.maxcount=100;
      param1.gamma=0.04;
      param1.epsilon=0.0001;
      param1.alpha=1.01;
      param1.verbose=false;

      
      // construct CGNEAlg
      ChebPolicy<T> p1;
      p1.assign(param1);
      est.zero();
      ChebAlg<T> * iter1 = p1.build(est,lop,rhs,rnorm,nrnorm,strout);
      iter1->run();
      
      // compute naive grad
      lop.applyOp(est,dltd);
      dltd.linComb(-1.0,rhs);
      sblop.applyAdjOp(est,dltd,grad);
      
      
      // construct LinFitLS function
      LinFitLS<T, ChebPolicy<T>, ChebPolicyData<T> > f1(nlop,B,rhs,param1,strout);
      
      FunctionalEvaluation<T> feval1(f1,x0);
      strout<<"\n Results from running Chebyshev\n";
      strout<<"value from LinFitLS    = "<<feval1.getValue()<<endl;
      strout<<"rnorm from Cheb        = "<<rnorm<<endl;
      strout<<"rnorm from LinFitLS    = "<<sqrt(2.0*abs(feval1.getValue()))<<endl;
      strout<<"difference rnorms      = "<<rnorm-sqrt(2.0*abs(feval1.getValue()))<<endl;
      // compute correction term of grad
      ChebPolicy<T> p2;
      p2.assign(param1);
      est.zero();
      ChebAlg<T> * iter2 = p2.build(est,lop,dltd,rnorm,nrnorm,strout);
      iter2->run();
      Vector<T> tmp(dom);
      sblop.applyAdjOp(est,rhs,tmp);
      grad.linComb(1.0, tmp);
      strout << "Output Gradient Norm:\n";
      strout << "Grad norm from LinFitLS         = "<< feval1.getGradient().norm() << endl;
      strout << "Grad norm computed in this test = "<<grad.norm() << endl;
      strout<<"*******************************************************"<<endl;

    strout.close();

    delete iter;
  }
  catch (RVLException & e) {
    e.write(str);
    exit(1);
  }
}

int main() {

  try {
    PlantSeeds(19490615);
    int rows=261;
    int cols=100;

    LinFitTestNL<float>(rows,cols,0.0,100,0.001,0.001,1000000.0,"2.1","float","./cgneflt.rpt",cout);
    LinFitTestNL<double>(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.2","double","./cgnedbl.rpt",cout);
 //   LinFitTestNL<complex<float> >(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.3","complex<float>","./cgnefltcplx.rpt",cout);
    LinFitTestNL<double>(rows,cols,0.2,100,0.001,0.001,1000000.0,"2.4","complex<double>","./cgnedblcplx.rpt",cout);
    LinFitTestNL<float>(rows,cols,0.0,100,0.001,0.001,1.0,"2.5","float","./cgneflt1.rpt",cout);
    LinFitTestNL<double>(rows,cols,0.2,100,0.001,0.001,1.0,"2.6","complex<double>","./cgnedblcplx1.rpt",cout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}

