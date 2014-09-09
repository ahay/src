#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
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
    IWaveOp op(*pars,stream);

    Vector<float> m(op.getDomain());

    AssignFilename mfn(valparse<std::string>(*pars,"csq"));
    Components<float> cm(m);
    cm[0].eval(mfn);
    OperatorEvaluation<float> opeval(op,m);

    std::string dnm = valparse<std::string>(*pars,"csq_d1_in","");
    if (dnm=="") {
      std::stringstream res;
      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      AdjointTest<float>(opeval.getDeriv(),rnd,res);
      if (retrieveGlobalRank()==0) {
	std::string outfile = valparse<std::string>(*pars,"outfile","");
	if (outfile.size()>0) {
	  ofstream outf(outfile.c_str());
	  outf<<res.str();
	  outf.close();
	}
	else {
	  cout<<res.str();
	}
      }
    }
    else {
      Vector<float> dmin(op.getDomain());
      Vector<float> d(op.getRange());
      Vector<float> dmout(op.getDomain());
      AssignFilename dminfn(dnm);
      Components<float> cdmin(dmin);
      cdmin[0].eval(dminfn);
      std::string dmoutnm = valparse<std::string>(*pars,"csq_d1_out","");
      if (dmoutnm.size()>0) {
	AssignFilename dmoutfn(dmoutnm);
	Components<float> cdmout(dmout);
	cdmout[0].eval(dmoutfn);
      }
      std::string dnm = valparse<std::string>(*pars,"data_in","");
      if (dnm.size()>0) {
	AssignFilename dfn(dnm);
	Components<float> cd(d);
	cd[0].eval(dfn);
      }
      opeval.getDeriv().applyOp(dmin,d);
      opeval.getDeriv().applyAdjOp(d,dmout);

      cerr<<"model space inner product = "<<dmin.inner(dmout)<<"\n";
      cerr<<"data space inner product  = "<<d.inner(d)<<"\n";
    }
    
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
