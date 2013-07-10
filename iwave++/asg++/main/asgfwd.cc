#include "asg_gfdm.h"

#include "asg_sampler.hh"
#include "seamx_headers.hh"
#include "iwop.hh"
#include "state.hh"
#include "samp.hh"
#include "sim.hh"
#include "pol.hh"

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

#include "parserpp.hh"

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
using TSOpt::StdSimData;
using TSOpt::OpNewCreatePolicy;
using TSOpt::PolicyBase;

using TSOpt::CPSim;

#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using TSOpt::GridWindowOp;
using TSOpt::SEGYLinMute;

using RVL::StdProductSpace;
using RVL::Vector;
using RVL::Components;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::Operator;
using RVL::OpComp;
using RVL::OpFO;
using RVL::LNLOperator;
using RVL::OperatorEvaluation;
using RVL::RVLException;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::parse;
using RVL::parse_except;

namespace ASG {
   typedef OpNewCreatePolicy<IWaveState,ASGSampler> ASGSamplerPolicy;
   typedef OpNewCreatePolicy<IWaveLinState,ASGLinSampler> ASGLinSamplerPolicy;

  typedef OpNewCreatePolicy<StdSimData<IWaveState>, CPSim<IWaveState,TSIndex> > FwdCPSimPolicy;

  /* dummy sampler policies to fill in the rest of the list */
  class LinSamplerPolicy: public PolicyBase<IWaveState,LinSampler > {};
  class AdjSamplerPolicy: public PolicyBase<IWaveState,LinSampler > {};

  typedef OpNewCreatePolicy<StdSimData<IWaveState>, StdSim<IWaveState> > StdIWavePolicy;
  
  /* dummy sim policies to fill in the rest of the list */
  class FwdRCSimPolicy: public PolicyBase<StdSimData<IWaveState>, StdSim<IWaveState> > {};
  class LinSimPolicy: public PolicyBase<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > {};
  class AdjFwdSimPolicy: public PolicyBase<StdSimData<IWaveState>, StdSim<IWaveState> > {};
  class AdjSimPolicy: public PolicyBase<StdSimData<IWaveLinState>, StdSim<IWaveLinState> > {};
} 

char ** xargv;

int main(int argc, char ** argv) {

  using namespace ASG;
  try { 

    // put there to avoid forgetting
    if (argc<2) { 
      cerr<<"asgfwd: fwd map\n";
      cerr<<"usage: asgfwd.x par=<par file>\n";
      exit(0);
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
      fflush(stream);
    }
    else {
#endif

      if (retrieveRank()==0) { cerr<<"IWAVE++::ASG forward map\n"; }

      /*********************************************************
       *             PARAMETER EXTRACTION                      *
       *********************************************************/

      string bulkname = "";
      string buoyname = "";
      string hdrname = "";
      string trcname = "";
      parse_except<string>(*pars,"bulkmod",bulkname);
      parse_except<string>(*pars,"buoyancy",buoyname);
      parse_except<string>(*pars,"datafile",trcname);
      parse_except<string>(*pars,"hdrfile",hdrname);

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
      MPIGridSpace m1sp(bulkname, "bulkmod");
      MPIGridSpace m2sp(buoyname, "buoyancy");
#else
      GridSpace m1sp(bulkname, "bulkmod");
      GridSpace m2sp(buoyname, "buoyancy");
#endif

      StdProductSpace<float> msp(m1sp,m2sp);
      Vector<float> x(msp);
      Components<float> cm(x);

      /* make SEGY space and vector */
#ifdef IWAVE_USE_MPI
      MPISEGYSpace tsp(hdrname);
#else
      SEGYSpace tsp(hdrname);
#endif
      Vector<float> y(tsp);

      /* assign files */
      AssignFilename m1fn(bulkname);
      AssignFilename m2fn(buoyname);
      AssignFilename tfn(trcname);

      cm[0].eval(m1fn);
      cm[1].eval(m2fn);
      y.eval(tfn);
	
      /*********************************************************
       *     OPERATOR CONSTRUCCTION, EVALUATION                *
       *********************************************************/

      IWaveOp<
      ASGSamplerPolicy,
	LinSamplerPolicy,
	AdjSamplerPolicy,
	StdIWavePolicy, // FwdCPSimPolicy, 
	FwdRCSimPolicy,
	LinSimPolicy,
	AdjFwdSimPolicy,
	AdjSimPolicy
	> 
	op_1(msp,tsp,*pars,stream,asg_gfdm);

      SEGYLinMute mute(sm,tm,wm);
      LinearOpFO<float> op_2(tsp,tsp,mute,mute);
      
      OpComp<ireal> op(op_1, op_2);

      OperatorEvaluation<ireal> opeval(op,x);
	
      y.copy(opeval.getValue());

      /*********************************************************
       *                    CLEANUP                            *
       *********************************************************/
       
      iwave_fdestroy();

#ifdef IWAVE_USE_MPI
    /* end nontriv comm branch */
    }
    /*
    fprintf(stream,"at MPI_Barrier: MPI_COMM_WORLD = %d\n",MPI_COMM_WORLD);
    fflush(stream);
    cerr<<"before rk = "<<retrieveGlobalRank()<<endl;
    */
    MPI_Barrier(MPI_COMM_WORLD);
    //    cerr<<"after rk = "<<retrieveGlobalRank()<<endl;
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e<<"asgfwd.x: ABORT\n";
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}


