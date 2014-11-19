#ifndef __TSOPT_HELMFFTW__OPS__
#define __TSOPT_HELMFFTW__OPS__

#include "rn.hh"
#include "op.hh"
#include "productspace.hh"
#include "mpiserialfo.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include <fftw3.h>

using RVL::ScalarFieldTraits;
using RVL::SpaceTest;
using RVL::Operator;
using RVL::LinearOp;
using RVL::Space;
using RVL::ProductSpace;
using RVL::Vector;
using RVL::Components;
using RVL::ProtectedDivision;
using RVL::RnArray;
using RVL::RVLScale;

namespace TSOpt {
    
    using RVL::BinaryLocalFunctionObject;
    using RVL::RVLException;
    using RVL::ContentPackage;
    using RVL::LocalDataContainer;
	
    
    
    class HelmFFTWFO: public BinaryLocalFunctionObject<ireal> {
        
    private:
        RPNT weights;
        ireal power, datum;
        IPNT n_arr;
        RPNT d_arr;
        IPNT sbc;      // starting boundary condition for each dimension
        IPNT ebc;      // end boundary condition for each dimension
        HelmFFTWFO();
        
    public:
        HelmFFTWFO(IPNT const & _narr,
                   RPNT const & _darr,
                   RPNT const & _weights,
                   IPNT const & _sbc,
                   IPNT const & _ebc,
                   ireal _power=0.0f,
                   ireal _datum=0.0f
                   )
        : power(_power), datum(_datum) {
            IASN(n_arr,_narr);
            RASN(d_arr,_darr);
            RASN(weights,_weights);
            IASN(sbc, _sbc);
            IASN(ebc, _ebc);
        }
        
        HelmFFTWFO(HelmFFTWFO const & f)
        : power(f.power), datum(f.datum) {
            IASN(n_arr,f.n_arr);
            RASN(d_arr,f.d_arr);
            RASN(weights,f.weights);
            IASN(sbc,f.sbc);
            IASN(ebc,f.ebc);
        }
        
        using RVL::LocalEvaluation<ireal>::operator();
        void operator()(LocalDataContainer<ireal> & x,
                        LocalDataContainer<ireal> const & y);
        
        string getName() const { string tmp = "HelmFFTWFO"; return tmp; }
        
    };
    
    class GridHelmFFTWOp: public LinearOp<float> {
    private:
        
        Space<float> const & dom;
        RPNT weights;
        float power, datum;
        IPNT sbc;      // starting boundary condition for each dimension
        IPNT ebc;      // end boundary condition for each dimension
        
        // default construction disabled
        GridHelmFFTWOp();
        
    protected:
        
        void apply(const Vector<float> & x,
                   Vector<float> & y) const;
        
        void applyAdj(const Vector<float> & x,
                      Vector<float> & y) const;
        
    public:
        
        GridHelmFFTWOp(GridHelmFFTWOp const & A)
        : dom(A.dom), power(A.power), datum(A.datum) {
            RASN(weights,A.weights);
            IASN(sbc,A.sbc);
            IASN(ebc,A.ebc);
        }
        
        GridHelmFFTWOp(Space<float> const & _dom,
                       RPNT _weights,
                       IPNT _sbc,
                       IPNT _ebc,
                       float _power=0.0f,
                       float _datum=0.0f):
        dom(_dom), power(_power), datum(_datum){
            try{
                RASN(weights,_weights);
                IASN(sbc,_sbc);
                IASN(ebc,_ebc);
            }
            catch (RVLException & e) {
                e<<"\ncalled from GridHelmFFTWOp constructor\n";
                throw e;
            }
        }
        
        ~GridHelmFFTWOp() {}
        
        // this class is considered terminal, with no overrides foreseen,
        // so clone method is not virtual
        LinearOp<float> * clone() const { return new GridHelmFFTWOp(*this); }
        
        // access to domain, range
        const Space<float> & getDomain() const { return dom; }
        const Space<float> & getRange() const { return dom; }
        
        ostream & write(ostream & str) const {
            str<<"GridHelmFFTWOp\n";
            return str;
        }
        
    };
}
#endif
