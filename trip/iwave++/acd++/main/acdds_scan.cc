#include "acd_gfdm.h"

#include "state.hh"
#include "acd_sampler.hh"
#include "seamx_headers.hh"
#include "iwop.hh"
#include "samp.hh"
#include "sim.hh"
#include "pol.hh"
#include "blockop.hh"
#include "parserpp.hh"
#include "cgnealg.hh"
#include "LinFitLS.hh"
#include "scantest.hh"
#include "CPsim.hh"

#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#include "mpigridpp.hh"
#else
#include "segypp.hh"
#include "gridpp.hh"
#endif
#include "segyops.hh"
#include "gridops.hh"

using namespace RVL;

using TSOpt::IWaveEnvironment;
using TSOpt::IWaveState;
using TSOpt::IWaveLinState;
using TSOpt::IWaveStep;
using TSOpt::IWaveStaticInit;
using TSOpt::IWaveDynamicInit;
using TSOpt::IWaveOp;
using TSOpt::Sampler;
using TSOpt::LinSampler;
using TSOpt::Sim;
using TSOpt::StdSim;
using TSOpt::StdRCSim;
using TSOpt::StdSimData;

using TSOpt::CPSim;

#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using TSOpt::SEGYLinMute;
using TSOpt::GridWindowOp;
using TSOpt::GridDerivOp;

using TSOpt::OpNewCreatePolicy;
using TSOpt::PolicyBase;
using RVL::StdProductSpace;
using RVL::Vector;
using RVL::Components;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::Operator;
using RVL::OpComp;
using RVL::OpFO;
using RVL::InjectOp;
using RVL::TensorOp;
using RVL::OperatorEvaluation;
using RVL::DerivEvaluation;
using RVL::RVLException;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::parse;
using RVL::Scan;
using RVLUmin::CGNEAlg;
using RVLUmin::CGNEPolicy;
using RVLUmin::CGNEPolicyData;
using RVLUmin::LinFitLS;

namespace ACD{

  /* dummy sampler policies to fill in the rest of the list */
  // class LinSamplerPolicy: public PolicyBase<IWaveState,LinSampler<IWaveLinState> > {};
  // class AdjSamplerPolicy: public PolicyBase<IWaveState,LinSampler<IWaveLinState> > {};
 
  /* Sampler Policies */
  /** FwdSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveState,ACDSampler> ACDSamplerPolicy;
  /** LinSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveLinState,ACDLinSampler> ACDLinSamplerPolicy;
  /** AdjSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveLinState,ACDAdjSampler> ACDAdjSamplerPolicy;
 
  /* Sim Policies */  
  /** FwdSimPolicy */
  typedef OpNewCreatePolicy<StdSimData<IWaveState>, StdSim<IWaveState> > StdIWavePolicy;
  /** LinFwdSimPolicy */
  typedef OpNewCreatePolicy<StdSimData<IWaveState>, StdRCSim<IWaveState> > StdRCIWavePolicy; 
  /** AdjFwdSimPolicy */
  typedef OpNewCreatePolicy<StdSimData<IWaveState>, CPSim<IWaveState,TSIndex> > FwdCPSimPolicy;
  /** LinSimPolicy and AdjSimPolicy */
  typedef OpNewCreatePolicy<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > LinSimPolicy;
 
  //typedef OpNewCreatePolicy<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > AdjSimPolicy;
   
  /* dummy sim policies to fill in the rest of the list */
  //class LinFwdSimPolicy: public PolicyBase<StdSimData<IWaveState>, StdSim<IWaveState> > {};
  //class LinSimPolicy: public PolicyBase<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > {};
  //class AdjFwdSimPolicy: public PolicyBase<StdSimData<IWaveState>, StdSim<IWaveState> > {};
  class AdjSimPolicy: public PolicyBase<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > {};
}

char ** xargv;

int main(int argc, char ** argv) {

  using namespace ACD;
  try { 

    // put there to avoid forgetting
    if (argc<2) { 
      RVLException e;
      e<<"acdadj: adjoint action of linearized fwd map\n";
      e<<"usage: acdadj.x par=<par file>\n";
      throw e;
    }

    /* set up execution environment */
    int ts=0;
#ifdef IWAVE_USE_MPI
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);
#endif    
    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc,argv,ts,&pars,&stream);

#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
    }
    else {
#endif

      /*********************************************************
       *             PARAMETER EXTRACTION                      *
       *********************************************************/

      /* files for acoustic params */
      string csqbegname = "";
      string csqendname = "";
      string dcsqname   = "";
      string hdrname = "";
      string trcname = "";
      string outfile = "";
      parse<string>(*pars,"outfile",outfile);
      parse_except<string>(*pars,"csq_beg",csqbegname);
      parse_except<string>(*pars,"csq_end",csqendname);
      parse_except<string>(*pars,"dcsq",dcsqname);
      parse_except<string>(*pars,"hdrfile",hdrname);
      parse_except<string>(*pars,"datafile",trcname);

      // assign window widths - default = 0;
      RPNT w;
      RASN(w,RPNT_0);
      parse<ireal>(*pars,"ww1",w[0]);
      parse<ireal>(*pars,"ww2",w[1]);
      parse<ireal>(*pars,"ww3",w[2]);

      ireal sm=0.0f;
      ireal wm=0.0f;
      ireal tm=0.0f;
      parse<ireal>(*pars,"mute_slope",sm);
      parse<ireal>(*pars,"mute_zotime",tm);
      parse<ireal>(*pars,"mute_width",wm);

      ireal rtol = 100.0*numeric_limits<ireal>::epsilon();
      ireal nrtol = 100.0*numeric_limits<ireal>::epsilon();
      int maxcount = 10;
      ireal maxstep = numeric_limits<ireal>::max();
      parse<ireal>(*pars,"ResidualTol",rtol);
      parse<ireal>(*pars,"GradientTol",nrtol);
      parse<int>(*pars,"MaxIter",maxcount);
      parse<ireal>(*pars,"MaxStep",maxstep);

      ireal alpha = REAL_ZERO;
      parse<ireal>(*pars,"DSParam",alpha);

      int nscan = 1;
      ireal hmin = REAL_ZERO;
      ireal hmax = REAL_ONE;
      int print_solver = 0;
      parse<int>(*pars,"NScan",nscan);
      parse<ireal>(*pars,"HMin",hmin);
      parse<ireal>(*pars,"HMax",hmax);
      parse<int>(*pars,"PrintSolver",print_solver);

      /*********************************************************
       *               INPUT CONSTRUCTION                      *
       *********************************************************/

      // construct spaces - note that pert space dmsp 
      // conveys window info
#ifdef IWAVE_USE_MPI
      MPIGridSpace msp(csqbegname, "csq", true);
      MPIGridSpace dmsp(dcsqname, "dcsq", true);
#else
      GridSpace msp(csqbegname, "csq", true);
      GridSpace dmsp(dcsqname, "dcsq", true);  // incore data
#endif

      // trace space
#ifdef IWAVE_USE_MPI
      MPISEGYSpace tsp(hdrname);
#else
      SEGYSpace tsp(hdrname);
#endif
    
      // start point
      Vector<ireal> x_beg(msp);
      AssignFilename mfnb(csqbegname);
      x_beg.eval(mfnb);
      Vector<ireal> x_end(msp);
      AssignFilename mfne(csqendname);
      x_end.eval(mfne);
      Vector<ireal> x_diff(msp);
      x_diff.copy(x_end);
      x_diff.linComb(-1.0f,x_beg);

      // scan direction - initialization below      
      Vector<ireal> dx(dmsp);

      /*********************************************************
       *     OPERATOR CONSTRUCCTION                            *
       *********************************************************/

      /* simulator */
      IWaveOp<
      ACDSamplerPolicy,
	ACDLinSamplerPolicy,
	ACDAdjSamplerPolicy,
	StdIWavePolicy,
	StdRCIWavePolicy,
	LinSimPolicy,
	FwdCPSimPolicy,
	AdjSimPolicy
	> 
	iwop(msp,tsp,*pars,stream,acd_gfdm);

      GridWindowOp wop(dmsp,x_beg,w);
      SEGYLinMute mute(sm,tm,wm);
      LinearOpFO<ireal> mop(tsp,tsp,mute,mute);
      OpComp<ireal> op(iwop,mop);
      int ddim = 0;
      if (retrieveGlobalRank() == 0) {
	if (msp.getGrid().gdim - msp.getGrid().dim < 1) {
	  RVLException e;
	  e<<"Error: acdds.scan.cc\n";
	  e<<"  velocity space does not have extra axis for DS\n";
	  throw e;
	}
	ddim=msp.getGrid().dim;
      }
      GridDerivOp ds(msp,ddim,alpha);
      TensorOp<ireal> dsop(op,ds);

      /*********************************************************
       *     SCAN                                              *
       *********************************************************/

      // least squares target vector - [data,zero]^T
      Vector<ireal> y(dsop.getRange());
      Components<ireal> cy(y);
      AssignFilename tfn(trcname);
      cy[0].eval(tfn);
      cy[1].zero();
      
      /* work vector - displaced model */
      Vector<ireal> x0(wop.getDomain());
      x0.zero();

      /* project direction into domain */
      OperatorEvaluation<ireal> wopeval(wop,x0);
      wopeval.getDeriv().applyAdjOp(x_diff,dx);
      
      /* output stream */
      std::stringstream res;
      res<<scientific;

      /* least squares fit function - CG version */
      CGNEPolicyData<ireal> pd(rtol,nrtol,maxstep,maxcount,print_solver);
      LinFitLS<ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > lsfit(dsop,wopeval.getDeriv(),y,pd,res);
      
      /* scan */
      Scan(lsfit,x_beg,x_diff,nscan,hmin,hmax,res);

      /* dump to file or stdout */
      if (retrieveRank() == 0) {
	if (outfile.size()>0) {
	  ofstream outf(outfile.c_str());
	  outf<<res.str();
	  outf.close();
	}
	else {
	  cout<<res.str();
	}
      }

      /*********************************************************
       *                       CLEANUP                         *
       *********************************************************/

      iwave_fdestroy();  

#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif  

  }
  catch (RVLException & e) {
    e<<"acdds_scan.x: ABORT\n";
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
