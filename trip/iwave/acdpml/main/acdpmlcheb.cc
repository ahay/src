#include "acdpml_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "chebalg.hh"
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
using TSOpt::SEGYLinMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::ChebAlg;

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
    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
                    valparse<float>(*pars,"mute_zotime",0.0f),
                    valparse<float>(*pars,"mute_width",0.0f));
      
    LinearOpFO<float> muteop(iwop.getRange(),iwop.getRange(),mute,mute);
    OpComp<float> op(iwop,muteop);

    Vector<ireal> m(op.getDomain());
    Vector<ireal> dm(op.getDomain());
    Vector<ireal> dd(op.getRange());

    AssignFilename mfn(valparse<std::string>(*pars,"rcsq"));
    Components<ireal> cm(m);
    cm[0].eval(mfn);

    AssignFilename dmfn(valparse<std::string>(*pars,"icsq"));
    Components<ireal> cdm(dm);
    cdm[0].eval(dmfn);
    dm.zero();

    AssignFilename ddfn(valparse<std::string>(*pars,"data"));
    Components<ireal> cdd(dd);
    cdd[0].eval(ddfn);

    int maxcount=valparse<int>(*pars,"MaxIter",10);
    float gamma=valparse<float>(*pars,"gamma",0.04f);
    float epsilon=valparse<float>(*pars,"epsilon",0.1);
    float alpha=valparse<float>(*pars,"alpha",1.01);

    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    
    /* output stream */
    std::stringstream res;
    res<<scientific;
    
    res<<endl<<"*******************************************************"<<endl;
    res<<"* Acoustic Constant Density Linearized Inversion via";
    res<<"* Chebyshev Algorithm for Normal Eqns"<<endl;
    res<<"* max iterations    = "<<maxcount<<endl;
    res<<"* error reduction   = "<<epsilon<<endl;
    res<<"* inversion level   = "<<gamma<<endl;
    res<<"* fudge factor      = "<<alpha<<endl;
    res<<"*******************************************************"<<endl;
    
    /* create Cheb object */
    float rnorm;
    float nrnorm;
    OperatorEvaluation<ireal> opeval(op,m);
    AdjointTest<float>(opeval.getDeriv(),rnd,cerr);
    ChebAlg<float> alg(dm,opeval.getDeriv(),dd,
		       rnorm, nrnorm, gamma, epsilon, alpha, maxcount, res);
    float nrnorm0=nrnorm;
    float rnorm0=rnorm;
    
    alg.run();
    
    // display results
    res<<"\n ******* summary ********  "<<endl;
    res<<"initial residual norm      = "<<rnorm0<<endl;
    res<<"residual norm              = "<<rnorm<<endl;
    res<<"residual redn              = "<<rnorm/rnorm0<<endl;
    res<<"initial gradient norm      = "<<nrnorm0<<endl;
    res<<"gradient norm              = "<<nrnorm<<endl;
    res<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
    
    std::string dataest = valparse<std::string>(*pars,"dataest","");
    std::string datares = valparse<std::string>(*pars,"datares","");
    if (dataest.size()>0) {
      Vector<float> est(op.getRange());
      AssignFilename estfn(dataest);
      est.eval(estfn);
      opeval.getDeriv().applyOp(dm,est);
      if (datares.size()>0) {
	Vector<float> res(op.getRange());
	AssignFilename resfn(datares);
	res.eval(resfn);
	res.copy(dd);
	res.linComb(-1.0f,est);
      }
    }

    if (retrieveRank() == 0) {
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
