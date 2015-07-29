#ifndef __RVLALG_LINFIT_Cheb_H
#define __RVLALG_LINFIT_Cheb_H

/** Given an Operator F, a Vector d in the range of op, and a initial vector dx=x0
 implements the function
 \f$$
 f(x) = \inf_{dx} \|DF(x)dx - d\|^2
 \f$$
 as an RVL::Functional. The linear least squares solver is specified by
 policy.
 */

#include "chebalg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {
    
    using namespace RVL;
    using namespace RVLAlg;
    
    template<typename Scalar>
    class LinFitLSCheb: public Functional<Scalar> {
        
        typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
        
    private:
        
        Operator<Scalar> const & op;        // operator
        LinearOp<Scalar> const & preop;     // preconditioner
        LinearOp<Scalar> const & helmop;    // smoothing op applied to gradient
        Vector<Scalar> const & d;           // data
        Vector<Scalar> const & x0;          // input initial linear solution
        ChebPolicyData<Scalar> s;
        mutable Scalar specbd;                      // spectrum bound of DF(x)
        bool refine;                        // refine as in Kern & Symes 1994
        mutable  Vector<Scalar> dx;         // preimage of linear solution
        mutable  Vector<Scalar> dltx;       // linear solution
        mutable bool applied;
        ostream & str;
        
    protected:
        
        void apply(const Vector<Scalar> & x,
                   Scalar & val) const {
            try {
                /*         if (applied) {
                 RVLException e;
                 e<<"Error: LinFitLS::apply(x,val)\n";
                 e<<"already applied, may not alter\n";
                 throw e;
                 }
                 */
                atype rnorm;
                atype nrnorm;
                // access Operator through OperatorEvaluation
                OperatorEvaluation<Scalar> opeval(op,x);
                
                // Get Derivative of Operator
                LinearOp<Scalar> const & lop = opeval.getDeriv();
                
                // Composition of lop and preop
                OpComp<Scalar> gop(preop,lop);
                
                Vector<Scalar> tmp(gop.getDomain());
                tmp.zero();
                
                Vector<Scalar> d0(lop.getRange());
                lop.applyOp(x0,d0);
                d0.linComb(1.0,d,-1.0);
                
                dx.zero();
                // build least square solver , solve for dx
                OperatorEvaluation<Scalar> gopeval(gop,x0);
                
                ChebPolicy<Scalar> chebp;
                chebp.assign(s);
                ChebAlg<Scalar> * solver
                = chebp.build(dx,gopeval.getDeriv(),d0,rnorm,nrnorm,str);
                
                solver->run();
                specbd = solver-> getSpectrumBound();
                // get the value of objective function
                val = 0.5*rnorm*rnorm;
                
                preop.applyOp(dx,dltx);
                dltx.linComb(1.0,x0);
                
                applied = true;
                delete solver;
            }
            catch (RVLException & e) {
                e<<"\ncalled from LinFitLSCheb::apply\n";
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
                OperatorEvaluation<Scalar> opeval(op,x);
                LinearOp<Scalar> const & lop = opeval.getDeriv();
                SymmetricBilinearOp<Scalar> const & sblop = opeval.getDeriv2();
                Vector<Scalar> dltd(lop.getRange());
                Vector<Scalar> gtmp(g.getSpace());
                // compute dltx and dltd = DF * dltx - d
                lop.applyOp(dltx,dltd);
                //  cerr << "\n dltx.norm()=" << dltx.norm() << endl;
                dltd.linComb(-1.0,d);
                //cerr << "\n dltd.norm()=" << dltd.norm() << endl;
                // naive computation of gradient
                sblop.applyAdjOp(dltx,dltd,gtmp);
                //  cerr << "\n g.norm()=" << g.norm() << endl;
                
                // compute and add correction term to gradient
                // cerr << "\n LinFitLS refine=" << refine << endl;
                if (refine) {
                    // cerr << "\n LinFitLS: inside refine branch\n" ;
                    atype rnorm;
                    atype nrnorm;
                    OpComp<Scalar> gop(preop,lop);
                    Vector<Scalar> tmp(gop.getDomain());
                    Vector<Scalar> dx1(gop.getDomain());
                    tmp.zero();
                    dx1.zero();
                    OperatorEvaluation<Scalar> gopeval(gop,x0);
                    // solve DF * dx = dltd in LS sense
                    ChebPolicyData<Scalar> param;
                    param.maxcount=s.maxcount;
                    param.gamma=s.gamma;
                    param.epsilon=s.epsilon;
                    param.alpha=s.alpha;
                    param.lbd_est=specbd;
                    cerr << "\n estimated spectrum bound is " << specbd <<  endl;
                    param.verbose=false;
                    ChebPolicy<Scalar> chebp;
                    chebp.assign(param);
                    ChebAlg<Scalar> * solver = chebp.build(dx1,gopeval.getDeriv(),dltd,rnorm,nrnorm,str);
                    solver->run();
                    delete solver;
                    
                    Vector<Scalar> tmp2(g.getSpace());
                    Vector<Scalar> dx2(preop.getRange());
                    preop.applyOp(dx1,dx2);
                    // compute and add correction term tmp to gradient g
                    sblop.applyAdjOp(dx2,d,tmp2);
                    // cerr << "\n LinFitLS:applyGradient  term1.norm = " << g.norm() << endl;
                    // cerr << "\n LinFitLS:applyGradient  term2.norm = " << tmp.norm() << endl;
                    gtmp.linComb(1.0, tmp2);
                    // cerr << "\n LinFitLS:applyGradient  g.norm = " << g.norm() << endl;
                }
                helmop.applyOp(gtmp,g);
            }
            catch (RVLException & e) {
                e<<"\ncalled from LinFitLS::applyGradient\n";
                throw e;
            }
            
        }
        
        void applyHessian(const Vector<Scalar> & x,
                          const Vector<Scalar> & dx,
                          Vector<Scalar> & dy) const {}
        
        Functional<Scalar> * clone() const {
            return new LinFitLSCheb<Scalar>(*this);
        }
        
    public:
        
        /* typical policy data 
         atype _rtol,
         atype _nrtol,
         int _maxcount,
         */
        LinFitLSCheb(Operator<Scalar> const & _op,
                     LinearOp<Scalar> const & _preop,
                     LinearOp<Scalar> const & _helmop,
                     Vector<Scalar> const & _d,
                     Vector<Scalar> const & _x0,
                     ChebPolicyData<Scalar> const & _s,
                     bool _refine=false,
                     ostream & _str=cerr)
        : op(_op),  helmop(_helmop), preop(_preop), d(_d), x0(_x0),s(_s),
        refine(_refine), dx(preop.getDomain()), dltx(preop.getRange()), 
        applied(false), str(_str) {
            try{
                dx.zero();
                if (s.verbose) {
                    str<<"\n";
                    str<<"==============================================\n";
                    str<<"LinFitLSCheb constructor - ls policy data = \n";
                    s.write(str);
                }
            }
            catch (RVLException & e) {
                e<<"\ncalled from LinFitLSCheb::Constructor\n";
                throw e;
            }
        }
        
        LinFitLSCheb(LinFitLSCheb<Scalar> const & f)
        : op(f.op), preop(f.preop), helmop(f.helmop), d(f.d), x0(f.x0), s(f.s), refine(f.refine),
        dx(f.dx), dltx(f.dltx), applied(f.applied), str(f.str) {}
        
        const Space<Scalar> & getDomain() const { return op.getDomain(); }
        
        Scalar getMaxStep(const Vector<Scalar> & x,
                          const Vector<Scalar> & dx) const {
            try {
                return op.getMaxStep(x,dx);
            }
            catch (RVLException & e) {
                e<<"\ncalled from LinFitLS::getMaxStep\n";
                throw e;
            }
        }
        
        Vector<Scalar> const & getLSSoln() const { return dltx; }
        
        ostream & write(ostream & str) const {
            str<<"LinFitLSCheb: \n";
            str<<"*** operator:\n";
            op.write(str);
            str<<"*** data vector:\n";
            d.write(str);
            return str;
        }
    };
}
#endif
