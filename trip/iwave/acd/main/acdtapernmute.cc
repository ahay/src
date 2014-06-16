#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "adjtest.hh"
#include "segyops.hh"

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
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYTaperMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

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
        
        // the Op
        IWaveOp iwop(*pars,stream);
        
        SEGYTaperMute tnm(valparse<float>(*pars,"mute_slope",0.0f),
                         valparse<float>(*pars,"mute_zotime",0.0f),
                         valparse<float>(*pars,"mute_width",0.0f),0,
                         valparse<float>(*pars,"min_gx",0.0f),
                         valparse<float>(*pars,"max_gx",numeric_limits<float>::max()),
                         valparse<float>(*pars,"taper_width",0.0f),0);
        
        LinearOpFO<float> tnmop(iwop.getRange(),iwop.getRange(),tnm,tnm);
        
        Vector<ireal> ddin(tnmop.getRange());
        Vector<ireal> ddout(tnmop.getRange());
        
        AssignFilename ddinfn(valparse<std::string>(*pars,"data"));
        Components<ireal> cddin(ddin);
        cddin[0].eval(ddinfn);
        
        AssignFilename ddoutfn(valparse<std::string>(*pars,"data_out"));
        Components<ireal> cddout(ddout);
        cddout[0].eval(ddoutfn);
        
        
        RVLRandomize<float> rnd(getpid(),-1.0,1.0);
        
        OperatorEvaluation<ireal> opeval(tnmop,ddin);
        AdjointTest<float>(opeval.getDeriv(),rnd,cerr);
        
        ddout.copy(opeval.getValue());
        
        
#ifdef IWAVE_USE_MPI
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
