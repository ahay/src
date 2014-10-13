#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "adjtest.hh"
//#include "productspace.hh"

#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#include "mpigridpp.hh"
#else
#include "segypp.hh"
#include "gridpp.hh"
#endif
#include "segyops.hh"
#include "gridops.hh"

#include "op.hh"
#include "blockop.hh"
#include <par.h>
#include "adjtest.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
    {"csqext",    0, true,  true },
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
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
using TSOpt::GridWindowOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif
using TSOpt::GridHelmOp;

using RVL::ScaleOpFwd;
using RVL::TensorOp;

using TSOpt::GridExtendOp;
using TSOpt::GridDerivOp;
using RVLUmin::CGNEAlg;

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
        // the Op
        IWaveOp iwop(*pars,stream);
        
        // assign window widths - default = 0;
        RPNT w_arr;
        RASN(w_arr,RPNT_1);
        w_arr[0]=valparse<float>(*pars,"scale1",1.0f);
        w_arr[1]=valparse<float>(*pars,"scale2",1.0f);
        ireal power=0.0f;
        ireal datum=0.0f;
        power=valparse<float>(*pars,"power",0.0f);
        datum=valparse<float>(*pars,"datum",0.0f);

        /* generate physical model space */
#ifdef IWAVE_USE_MPI
        MPIGridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#else
        GridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#endif
        // make it a product, so it's compatible with domain of op
        StdProductSpace<ireal> dom(csqsp);
        
        Vector<ireal> m_in(dom);
        Vector<ireal> m_out(dom);
        
        AssignFilename minfn(valparse<std::string>(*pars,"csqin"));
        Components<ireal> cmin(m_in);
        cmin[0].eval(minfn);
        
        AssignFilename moutfn(valparse<std::string>(*pars,"csqout"));
        Components<ireal> cmout(m_out);
        cmout[0].eval(moutfn);
        
        GridHelmOp hop(dom,w_arr,power,datum);
        
        hop.applyOp(m_in,m_out);
        
        
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
