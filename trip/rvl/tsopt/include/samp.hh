#ifndef __TSOPT_SAMP
#define __TSOPT_SAMP

#include "sim.hh"
#include "alg.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;

  class Sampler: public Writeable {
    
  public:
    
    Sampler() {}
    Sampler(Sampler const &) {}
    virtual ~Sampler() {}

    /** initialize internal data before each time loop */
    virtual bool init() = 0;
    /** algorithm to run before each time step */
    virtual Algorithm & pre()  const = 0;
    /** algorithm to run after each time step */
    virtual Algorithm & post() const = 0;
    /** terminator to control time loop */
    virtual TimeTerm & term() const = 0;
    /** internal to external after each time loop */
    virtual void flush() const = 0;
    /** return start time of time loop as Time object */
    virtual Time const & start() const = 0;
    /** return final time of time loop as Time object */
    virtual Time const & end() const = 0;
    /** returns number of records in output data - for multisim */
    virtual int getNumRecs() const = 0;
    /** returns current record number */
    virtual int getRecIndx() const = 0;

  };

  /** additional algs handle single set of tangent fields -
      accomodates both lin and adj simulations */
  class LinSampler: public Sampler {
    
  public:
    
    LinSampler() {}
    LinSampler(LinSampler const &) {}
    virtual ~LinSampler() {}

    /** algorithm to run before each time step */
    virtual Algorithm & linpre()  const = 0;
    /** algorithm to run after each time step */
    virtual Algorithm & linpost() const = 0;
	/** terminator to control lin-sim time loop */
    virtual TimeTerm & linterm() const = 0;
  };

  /*
  template<typename State>
  class SamplerPolicy {
  public:
    virtual Sampler<State> * create(State &) const = 0;
  };
  */

}

#endif
