#ifndef __RVLALG_MULTIFIT_L2_H
#define __RVLALG_MULTIFIT_L2_H

/** Given an Operator F and a Vector d in the range of op,
 implements the function
 \f$$
 f(x,dx) = \|DF(x)dx - d\|^2
 \f$$
 as an RVL::Functional. The linear least squares solver is specified by
 policy.
 */

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {
    
    using namespace RVL;
    using namespace RVLAlg;
    
    template<typename Scalar>
    class MultiFit: public Functional<Scalar> {
        
        typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
        
    private:
        
        StdProductSpace<Scalar> const & pdom;     // domain of this functional
        Operator<Scalar> const & op;        // operator
        LinearOp<Scalar> const & gextop;     // grid extension linear op
        LinearOp<Scalar> const & preop;     // preconditioner
        LinearOp<Scalar> const & helmop;    // smoothing op applied to gradient of background velocity
        Vector<Scalar> const & d;           // data
        mutable Vector<Scalar> dltd;        //
        mutable Vector<Scalar> extm;        // extended version of velocity
        mutable bool applied;
        ostream & str;
        
    protected:
        
        void apply(const Vector<Scalar> & x,
                   Scalar & val) const {
            try {
                Components<Scalar> cx(x);
                if (cx.getSize()!=2) {
                    RVLException e;
                    e << "Error: MultiFit::apply\n";
                    e<<"input data do not have two components \n";
                    throw e;
                }
                
                gextop.applyOp(cx[0],extm);

                // access Operator through OperatorEvaluation
                OperatorEvaluation<Scalar> opeval(op,extm);
                
                // Get Derivative of Operator
                LinearOp<Scalar> const & lop = opeval.getDeriv();
                
                // Composition of lop and preop
                OpComp<Scalar> gop(preop,lop);
                Vector<Scalar> tmp(gop.getDomain());
                tmp.zero();
                OperatorEvaluation<Scalar> gopeval(gop,tmp);
                
                gopeval.getDeriv().applyOp(cx[1],dltd);
                dltd.linComb(-1.0, d);
                // get the value of objective function
                val=dltd.normsq()/2.0;
                applied = true;
            }
            catch (RVLException & e) {
                e<<"\ncalled from MultiFit::apply\n";
                throw e;
            }
        }
        
        void applyGradient(const Vector<Scalar> & x,
                           Vector<Scalar> & g) const {
            try{
                if(!applied){
                    Scalar val;
                    this->apply(x,val);
                    // cerr << "\n val=" << val << endl;
                }
                
                Components<Scalar> cx(x);
                Components<Scalar> cg(g);
                if (cx.getSize()!=cg.getSize()) {
                    RVLException e;
                    e << "Error: MultiFit::apply\n";
                    e<<"model and gradient do not have the same dimension\n";
                    throw e;
                }
                //cerr << "\n extm.norm() = "<< extm.norm() << endl; 
                OperatorEvaluation<Scalar> opeval(op,extm);
                LinearOp<Scalar> const & lop = opeval.getDeriv();
                // Composition of lop and preop
                OpComp<Scalar> gop(preop,lop);
                Vector<Scalar> tmp(gop.getDomain());
                tmp.zero();
                OperatorEvaluation<Scalar> gopeval(gop,tmp);
                
                SymmetricBilinearOp<Scalar> const & sblop = opeval.getDeriv2();
                Vector<Scalar> tmp1(op.getDomain());
                tmp1.zero();
                Vector<Scalar> gtmp(cg[0].getSpace());
                str << "\n cx[1] = " << cx[1].norm() << endl;
                str << "\n dltd = " << dltd.norm() << endl;
                // computation of gradient of velocity
                sblop.applyAdjOp(cx[1],dltd,tmp1);
                str << "\n tmp1 = " << tmp1.norm() << endl;
                gextop.applyAdjOp(tmp1,gtmp);
                // compute gradient of reflectivity
                gopeval.getDeriv().applyAdjOp(dltd,cg[1]);
                str << "\n gtmp = " << gtmp.norm() << endl;
                str << "\n before helmop cg[0] = " << cg[0].norm() << endl;
                helmop.applyOp(gtmp,cg[0]);

                str << "\n after helmop cg[0] = " << cg[0].norm() << endl;
            }
            catch (RVLException & e) {
                e<<"\ncalled from MultiFit::applyGradient\n";
                throw e;
            }
            
        }
        
        void applyHessian(const Vector<Scalar> & x,
                          const Vector<Scalar> & dx,
                          Vector<Scalar> & dy) const {}
        
        Functional<Scalar> * clone() const {
            return new MultiFit<Scalar>(*this);
        }
        
    public:
        
        /* typical policy data
         atype _rtol,
         atype _nrtol,
         int _maxcount,
         */
        MultiFit(StdProductSpace<Scalar> const & _pdom,
                   Operator<Scalar> const & _op,
                   LinearOp<Scalar> const & _gextop,
                   LinearOp<Scalar> const & _preop,
                   LinearOp<Scalar> const & _helmop,
                   Vector<Scalar> const & _d,
                   ostream & _str=cerr)
        : pdom(_pdom), op(_op), gextop(_gextop), preop(_preop), helmop(_helmop), d(_d), dltd(op.getRange()),
        extm(op.getDomain()), applied(false), str(_str) {
            try{
                dltd.zero();
                extm.zero();
            }
            catch (RVLException & e) {
                e<<"\ncalled from MultiFit::Constructor\n";
                throw e;
            }
        }
        
        MultiFit(MultiFit<Scalar> const & f)
        : pdom(f.pdom), op(f.op), gextop(f.gextop), preop(f.preop), helmop(f.helmop), d(f.d),
        dltd(f.dltd), extm(f.extm), applied(f.applied), str(f.str) {}
        
        const Space<Scalar> & getDomain() const { return pdom; }
        
        Scalar getMaxStep(const Vector<Scalar> & x,
                          const Vector<Scalar> & dx) const {
            try {
                return op.getMaxStep(x,dx);
            }
            catch (RVLException & e) {
                e<<"\ncalled from MultiFit::getMaxStep\n";
                throw e;
            }
        }
        
        ostream & write(ostream & str) const {
            str<<"MultiFit: \n";
            str<<"*** operator:\n";
            op.write(str);
            str<<"*** data vector:\n";
            d.write(str);
            return str;
        }
    };
}
#endif
