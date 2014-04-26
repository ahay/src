    /* parameters for LBFGS:
	parameters:
	@param f - function to be minimized (RVL::Functional)
	@param x - solution RVL::Vector - initial guess on call, estimated solution on return
	@param _ihs - inverse Hessian scale - overall scale factor, so initial Hessian is this Scalar multiple of identity operator
	@param _maxits - max number of LBFGS iterations
	@param _mud - max stored BFGS updates - stored inverse Hessian approximation has this rank (at most)
	@param _maxsamp - max number of steps permitted in each line search
	@param _disp - verbosity flag - false = no output, true = function value, gradient norm at each iteration, report of line search
	@param _sl1 - first line search step
	@param _eta1 - lower G-A parameter
	@param _eta2 - upper G-A parameter
	@param _gamma1 - line search backtrack factor
	@param _gamma2 - line search extrapolation factor ("internal doubling")
	@param _maxfrac - fraction of max step to boundary permitted
	@param _agradtol - stopping tolerance for gradient norm, absolute
	@param _rgradtol - stopping tolerance for gradient norm, relative to initial gradient norm
	@param _str - verbose output unit
    */

#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "ls.hh"
#include "LBFGSBT.hh"
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
using RVL::valparse;
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
using RVL::RVLAssignConst;
using RVL::RVLRandomize;
using RVL::RVLMin;
using RVL::ULBoundsTest;
using RVL::FunctionalBd;
using RVL::StdLeastSquaresFcnlGN;
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
using RVL::MPISerialFunctionObjectRedn;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::LBFGSBT;

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
    
    IWaveOp iwop(*pars,stream);
      
    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
                       valparse<float>(*pars,"mute_zotime",0.0f),
                       valparse<float>(*pars,"mute_width",0.0f));
      
    LinearOpFO<float> muteop(iwop.getRange(),iwop.getRange(),mute,mute);
    OpComp<float> op(iwop,muteop);
    
    Vector<ireal> m(op.getDomain());
    Vector<ireal> dm(op.getDomain());
    Vector<ireal> dd(op.getRange());

    AssignFilename mfn(valparse<std::string>(*pars,"csq"));
    Components<ireal> cm(m);
    cm[0].eval(mfn);

    AssignFilename dmfn(valparse<std::string>(*pars,"icsq"));
    Components<ireal> cdm(dm);
    cdm[0].eval(dmfn);
    dm.copy(m);

    AssignFilename ddfn(valparse<std::string>(*pars,"data"));
    Components<ireal> cdd(dd);
    cdd[0].eval(ddfn);

    /* output stream */
    std::stringstream res;
    res<<scientific;
    
    StdLeastSquaresFcnlGN<float> fls(op,dd);

    // lower, upper bds for csq
    Vector<float> lb(op.getDomain());
    Vector<float> ub(op.getDomain());
    float lbcsq=valparse<float>(*pars,"cmin");
    lbcsq=lbcsq*lbcsq;
    RVLAssignConst<float> asl(lbcsq);
    lb.eval(asl);
    float ubcsq=valparse<float>(*pars,"cmax");
    ubcsq=ubcsq*ubcsq;
    RVLAssignConst<float> asu(ubcsq);
    ub.eval(asu);
    RVLMin<float> mn;
#ifdef IWAVE_USE_MPI
    MPISerialFunctionObjectRedn<float,float> mpimn(mn);
    ULBoundsTest<float> ultest(lb,ub,mpimn);
#else
    ULBoundsTest<float> ultest(lb,ub,mn);
#endif
    
    FunctionalBd<float> fbd(fls,ultest);

    LBFGSBT<float> alg(fbd,dm,
		       valparse<float>(*pars,"InvHessianScale",1.0f),
		       valparse<int>(*pars,"MaxInvHessianUpdates",5),
		       valparse<int>(*pars,"MaxLineSrchSteps",10),
		       valparse<bool>(*pars,"VerboseDisplay",true), 
		       valparse<float>(*pars,"FirstStepLength",1.0f), 
		       valparse<float>(*pars,"GAStepAcceptThresh",0.1f),
		       valparse<float>(*pars,"GAStepDoubleThresh",0.9f), 
		       valparse<float>(*pars,"LSBackTrackFac",0.5f), 
		       valparse<float>(*pars,"LSDoubleFac",1.8f),
		       valparse<float>(*pars,"MaxFracDistToBdry",1.0), 
		       valparse<float>(*pars,"LSMinStepFrac",1.e-06),
		       valparse<int>(*pars,"MaxLBFGSIter",3), 
		       valparse<float>(*pars,"AbsGradThresh",0.0), 
		       valparse<float>(*pars,"RelGradThresh",1.e-2), 
		       res);
    
    alg.run();
    
    if (retrieveRank() == 0) {
      std::string outfile = valparse<std::string>(*pars,"VerboseReport","");
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
