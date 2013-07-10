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
#include "adjtest.hh"

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
using RVL::OperatorEvaluation;
using RVL::DerivEvaluation;
using RVL::RVLException;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::AdjointTest;
using RVL::RVLRandomize;
using RVL::parse;
using RVL::parse_except;

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

      if (retrieveRank()==0) { cerr<<"IWAVE++::ACD adjoint test\n"; }

      /* files for acoustic params */
      string csqname = "";
      string dcsqname = "";
      string outfile = "";
      string hdrname = "";
      parse<string>(*pars,"outfile",outfile);
      parse_except<string>(*pars,"csq",csqname);
      parse_except<string>(*pars,"dcsq",dcsqname);
      parse_except<string>(*pars,"hdrfile",hdrname);

      // assign window widths - default = 0;
      RPNT w;
      RASN(w,RPNT_0);
      parse<float>(*pars,"ww1",w[0]);
      parse<float>(*pars,"ww2",w[1]);
      parse<float>(*pars,"ww3",w[2]);

      float sm=0.0f;
      float wm=0.0f;
      float tm=0.0f;
      parse<float>(*pars,"mute_slope",sm);
      parse<float>(*pars,"mute_zotime",tm);
      parse<float>(*pars,"mute_width",wm);

      /*********************************************************
       *               INPUT CONSTRUCTION                      *
       *********************************************************/

#ifdef IWAVE_USE_MPI
      MPIGridSpace msp(csqname, "csq");
      MPIGridSpace dmsp(dcsqname, "dcsq");
#else
      GridSpace msp(csqname, "csq");
      GridSpace dmsp(dcsqname, "dcsq");
#endif
      Vector<float> x(msp);

      /* make SEGY space and vector */
#ifdef IWAVE_USE_MPI
      MPISEGYSpace tsp(hdrname);
#else
      SEGYSpace tsp(hdrname);
#endif

      /* assign files */
      AssignFilename mfn(csqname);
      x.eval(mfn);

      /*********************************************************
       *     OPERATOR CONSTRUCCTION, EVALUATION                *
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

      GridWindowOp wop(dmsp,x,w);
      SEGYLinMute mute(sm,tm,wm);
      LinearOpFO<float> mop(tsp,tsp,mute,mute);
      OpComp<float> op(wop,iwop,mop);

      /* work vector - displaced model */
      Vector<float> x0(op.getDomain());
      x0.zero();

      /* evaluate at (displaced) origin */
      OperatorEvaluation<float> opeval(op,x0);

      RVLRandomize<float> rnd(getpid(),-1.0f, 1.0f);

      std::stringstream res;
      AdjointTest(opeval.getDeriv(),rnd,res);
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
    e<<"acdadj.x: ABORT\n";
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
