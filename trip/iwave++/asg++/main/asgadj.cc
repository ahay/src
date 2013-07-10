#include "asg_gfdm.h"

#include "state.hh"
#include "asg_sampler.hh"
#include "seamx_headers.hh"
#include "iwop.hh"
#include "samp.hh"
#include "sim.hh"
#include "pol.hh"
#include "blockop.hh"
#include "parserpp.hh"

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
using RVL::parse;

namespace ASG{

  /* dummy sampler policies to fill in the rest of the list */
  // class LinSamplerPolicy: public PolicyBase<IWaveState,LinSampler<IWaveLinState> > {};
  // class AdjSamplerPolicy: public PolicyBase<IWaveState,LinSampler<IWaveLinState> > {};
 
  /* Sampler Policies */
  /** FwdSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveState,ASGSampler> ASGSamplerPolicy;
  /** LinSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveLinState,ASGLinSampler> ASGLinSamplerPolicy;
  /** AdjSamplerPolicy */
  typedef OpNewCreatePolicy<IWaveLinState,ASGAdjSampler> ASGAdjSamplerPolicy;
 
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

  using namespace ASG;
  try { 

    // put there to avoid forgetting
    if (argc<2) { 
      RVLException e;
      e<<"asgadj: adjoint action of linearized fwd map\n";
      e<<"usage: asgadj.x par=<par file>\n";
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

      /* concoct acoustic model space, vector, product vector components */
      string bulkname = "";
      string buoyname = "";
      string dbulkname = "";
      string dbuoyname = "";
      string mbulkname = "";
      string mbuoyname = "";
      string hdrname = "";
      string trcname = "";
      parse_except<string>(*pars,"bulkmod",bulkname);
      parse_except<string>(*pars,"buoyancy",buoyname);
      parse_except<string>(*pars,"dbulkmod",dbulkname);
      parse_except<string>(*pars,"dbuoyancy",dbuoyname);
      parse_except<string>(*pars,"mbulkmod",mbulkname);
      parse_except<string>(*pars,"mbuoyancy",mbuoyname);
      parse_except<string>(*pars,"hdrfile",hdrname);
      parse_except<string>(*pars,"datafile",trcname);

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
      MPIGridSpace m1sp(bulkname, "bulkmod");
      MPIGridSpace dm1sp(dbulkname, "bulkmod");
      MPIGridSpace m2sp(buoyname, "buoyancy");
      MPIGridSpace dm2sp(dbuoyname, "buoyancy");
#else
      GridSpace m1sp(bulkname, "bulkmod");
      GridSpace dm1sp(dbulkname, "bulkmod");
      GridSpace m2sp(buoyname, "buoyancy");
      GridSpace dm2sp(dbuoyname, "buoyancy");
#endif
      StdProductSpace<float> msp(m1sp,m2sp);
      Vector<float> x(msp);
      Components<float> cx(x);

      StdProductSpace<float> dmsp(dm1sp,dm2sp);
      Vector<float> dx(dmsp);
      Components<float> cdx(dx);

      /* make SEGY space and vector */
#ifdef IWAVE_USE_MPI
      MPISEGYSpace tsp(hdrname);
#else
      SEGYSpace tsp(hdrname);
#endif
    
      Vector<float> dy(tsp);

      /* assign files */
      AssignFilename m1fn(bulkname);
      AssignFilename m2fn(buoyname);
      cx[0].eval(m1fn);
      cx[1].eval(m2fn);

      /* assign files */
      AssignFilename mm1fn(mbulkname);
      AssignFilename mm2fn(mbuoyname);
      cdx[0].eval(mm1fn);
      cdx[1].eval(mm2fn);

      AssignFilename tfn(trcname);
      dy.eval(tfn);

      /*********************************************************
       *     OPERATOR CONSTRUCCTION, EVALUATION                *
       *********************************************************/
      
      /* simulator */
      IWaveOp<
      ASGSamplerPolicy,
	ASGLinSamplerPolicy,
	ASGAdjSamplerPolicy,
	StdIWavePolicy,
	StdRCIWavePolicy,
	LinSimPolicy,
	FwdCPSimPolicy,
	AdjSimPolicy
	> 
	iwop(msp,tsp,*pars,stream,asg_gfdm);

      GridWindowOp wop(dmsp,x,w);

      SEGYLinMute mute(sm,tm,wm);
      LinearOpFO<float> mop(tsp,tsp,mute,mute);

      OpComp<float> op(wop, iwop, mop);

      Vector<float> x0(op.getDomain());
      x0.zero();

      /* evaluate at (displaced) origin */
      OperatorEvaluation<float> opeval(op,x0);

      /* evaluate derivative */    
      opeval.getDeriv().applyAdjOp(dy,dx);

      iwave_fdestroy();  

#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif  

  }
  catch (RVLException & e) {
    e<<"asgadj.x: ABORT\n";
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
