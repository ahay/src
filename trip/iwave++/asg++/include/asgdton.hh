#ifndef __ASG_DtoN__
#define __ASG_DtoN__

#include "asg_headers.hh"
#include "asg_sampler.hh"
#include "seamx_headers.hh"
#include "iwop.hh"
#include "state.hh"
#include "samp.hh"
#include "sim.hh"

#include "CPsim.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#else
#include "gridpp.hh"
#include "segypp.hh"
#endif

namespace ASG {

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
  using TSOpt::segytrace;

  using RVL::LocalDataContainer;
  using RVL::StdProductSpace;
  using RVL::Space;
  using RVL::Vector;
  using RVL::Space;
  using RVL::Components;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::OperatorEvaluation;
  using RVL::RVLException;
  using RVL::AssignFilename;
  using RVL::AssignTag;
  using RVL::AssignParams;
  using RVL::BinaryLocalFunctionObject;

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

  typedef IWaveOp<
      ASGSamplerPolicy,
      ASGLinSamplerPolicy,
      ASGAdjSamplerPolicy,
      StdIWavePolicy,
      StdRCIWavePolicy,
      LinSimPolicy,
      FwdCPSimPolicy,
      AdjSimPolicy
      > ASGSimOp;

  class FwdInt: public RVL::BinaryLocalFunctionObject<float> {
  public:
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in);

    string getName() const { string ret="FwdInt"; return ret; }
  };

  class AdjInt: public RVL::BinaryLocalFunctionObject<float> {
  public:
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in);

    string getName() const { string ret="AdjInt"; return ret; }
  };

  class TimeRev: public RVL::BinaryLocalFunctionObject<float> {
  public:
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in);

    string getName() const { string ret="TimeRev"; return ret; }
  };
    
  class ASGDtoN: public RVL::LinearOp<float> {

  private:

    FILE * str;

#ifdef IWAVE_USE_MPI
    MPISEGYSpace const & dom;
#else
    SEGYSpace const & dom;
#endif

    mutable Vector<float> w;     // domain workspace - trace intgl
    Vector<float> const & mdl;
    bool integrate;              // flag for trace integration
    
    ASGSimOp op;                 // underlying asg op

    ASGDtoN();

  protected:

    LinearOp<float> * clone() const {
      return new ASGDtoN(*this);
    }

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

  public:

    ASGDtoN(Vector<float> const & _mdl,
#ifdef IWAVE_USE_MPI
	    MPISEGYSpace const & _dom,
#else
	    SEGYSpace const & _dom,
#endif	  
	    PARARRAY const & _par,
	    FILE * _str);

    ASGDtoN(ASGDtoN const &);
    ~ASGDtoN() {}

    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    Operator<float> & getOp() { return op; }
    Operator<float> const & getOp() const { return op; }

    ostream & write(ostream & str) const;
  };

}

#endif
