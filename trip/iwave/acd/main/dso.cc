#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "gridops.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "adjtest.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
    {"csq",    0, true,  true },
    {"data",   1, false, true },
    {"source", 1, true,  false},
    {"",       0, false, false}
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::AdjointTest;
using TSOpt::IWaveEnvironment;
using TSOpt::GridWindowOp;
using TSOpt::GridMaskOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
#else
using TSOpt::GridSpace;
#endif
using TSOpt::GridDerivOp;


int xargc;
char **xargv;

int main(int argc, char ** argv) {
    
    try {
        
#ifdef IWAVE_USE_MPI
        int ts=0;
        MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);
#endif
        
        PARARRAY * pars = NULL;
        FILE * stream = NULL;
        IWaveEnvironment(argc, argv, 0, &pars, &stream);
#ifdef IWAVE_USE_MPI
        if (retrieveGroupID() == MPI_UNDEFINED) {
            fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
            fflush(stream);
        }
        else {
#endif
            
            
            /* generate physical model space */
#ifdef IWAVE_USE_MPI
            MPIGridSpace csqsp(valparse<std::string>(*pars,"csqext"),"notype",true);
#else
            GridSpace csqsp(valparse<std::string>(*pars,"csqext"),"notype",true);
#endif
            
            // make it a product, so it's compatible with domain of op
            StdProductSpace<ireal> dom(csqsp);
            int dsdir = INT_MAX;
            if (retrieveGlobalRank()==0) dsdir=csqsp.getGrid().dim;
#ifdef IWAVE_USE_MPI
            if (MPI_Bcast(&dsdir,1,MPI_INT,0,retrieveGlobalComm())) {
                RVLException e;
                e<<"Error: acddscheb_grad, rank="<<retrieveGlobalRank()<<"\n";
                e<<"  failed to bcast dsdir\n";
                throw e;
            }
#endif
            
            // assign window widths - default = 0;
            RPNT swind,ewind;
            RASN(swind,RPNT_0);
            RASN(ewind,RPNT_0);
            swind[0]=valparse<float>(*pars,"sww1",0.0f);
            swind[1]=valparse<float>(*pars,"sww2",0.0f);
            swind[2]=valparse<float>(*pars,"sww3",0.0f);
            ewind[0]=valparse<float>(*pars,"eww1",0.0f);
            ewind[1]=valparse<float>(*pars,"eww2",0.0f);
            ewind[2]=valparse<float>(*pars,"eww3",0.0f);
//            RASN(wind,RPNT_0);
//            wind[0]=valparse<float>(*pars,"ww1",0.0f);
//            wind[1]=valparse<float>(*pars,"ww2",0.0f);
//            wind[2]=valparse<float>(*pars,"ww3",0.0f);
            
            GridDerivOp dsop(dom,dsdir,valparse<float>(*pars,"DSWeight",1.0f));
            
            Vector<ireal> m_in(dom);
            Vector<ireal> tmp(dom);
            Vector<ireal> m_out(dom);
            
            AssignFilename m_infn(valparse<std::string>(*pars,"csqin"));
            Components<ireal> cm_in(m_in);
            cm_in[0].eval(m_infn);
            
            GridMaskOp wop(dom,m_in,swind,ewind);
            OperatorEvaluation<ireal> wopeval(wop,m_in);
            CompLinearOp<float> lop(wopeval.getDeriv(),dsop);
            
            AssignFilename m_outfn(valparse<std::string>(*pars,"csqout"));
            Components<ireal> cm_out(m_out);
            cm_out[0].eval(m_outfn);
            
            RVLRandomize<float> rnd(getpid(),-1.0,1.0);
            
            AdjointTest<float>(lop,rnd,cerr);
           
            tmp.zero();
            lop.applyOp(m_in,tmp);
            lop.applyOp(tmp,m_out);
            
#ifdef IWAVE_USE_MPI
        }
        MPI_Finalize();
#endif
    }
    catch (RVLException & e) {
        e.write(cerr);
#ifdef IWAVE_USE_MPI
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(1);
    }
    
}
