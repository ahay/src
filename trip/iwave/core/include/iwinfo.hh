#ifndef __IWAVE_M_TUPLE
#define __IWAVE_M_TUPLE

#include "std_cpp_includes.hh"
#include "utils.h"
#include "except.hh"
#include "write.hh"
#include "model.h"

/** This file defines the function interfaces for which
    implementations must be provided to create an IWAVE application.
 */ 

using RVL::RVLException;

size_t pow2(int);

/* max number of data tuples permitted */
#define IWAVEMAXDATA 256

// forward declaration
class IWaveInfo;

/**

Time step function. Interface for time steps of all orders of
derivatives - order inferred from structure of input. The first
argument is an array of \ref RDOM s - each \ref RDOM contains all
static and dynamic fields interacting in the simulation. For the kth
derivative (or adjoint derivative), 2<sup>k</sup> \ref RDOM s (copies
of the basic field setup) are required. The fwd flag indicates whether
the computation is part of a forward (true) or adjoint (false) time
loop. For multistep methods (for expl leapfrog), iv indicates the
substep of the overall step being computed.  <p>

An example design: this function can determine the order of derivative
by taking the log base 2 of dom.size(), extract the \ref RARR s from
the \ref RDOM s, the C arrays and array metadata from the \ref RARR
s, and pass the necessary information for time stepping to low-level
C functions.

 <p>
 @param[out] dom - std::vector<RDOM *> storing dynamic and static fields
 @param[in] fwd - flag for forward vs. backward (adjoint) time step
 @param[in] iv - internal (sub)step index
 @param[in] fdpars - opaque FD_MODEL parameter object 
 @return 0 on success, else error code

 called in IWaveSim::run

 */
typedef void (*FD_TIMESTEP)(std::vector<RDOM *> dom, bool fwd, int iv, void* fdpars);

/** FD model internals initializer - creates appropriate data
    structure to serve as fdpars parameter struct for \ref IMODEL
    specialization, initialize data members of \ref IMODEL.fdpars
    using data from pars, grid, and dt. The data members include (for
    example) scheme-characteristic finite difference coefficients,
    parameters and arrays for absorbing boundary conditions, arrays to
    indicate whether boundar faces are inter-domain or external for
    domain decomposition, and other parameters convenient for defining
    time steps.
  
    This function acts as the constructor for a concrete subclass of IMODEL.

    @param[in] pars - parameter array, assumed initialized.
    
    @param[in] stream - verbose output stream
    
    @param[in] g - primal simulation \ref grid , initialized via I/O
    on first IOKEY
    
    @param[in] dt - time step, will have been initialized on call by prior
    call to FD_TIMEGRID

    @param[out] specs - fd parameter struct containing app-particular
    info such as coefficient arrays, Courant numbers, PML damping
    arrays, etc. 
    
    @return - 0 on success, else error code.
*/

/*
typedef int (*FD_MODELINIT)(PARARRAY * pars, 
			    FILE * stream, 
			    grid const & g,
			    ireal dt,
			    std::vector<std::string> & active,
			    void ** specs);
*/
typedef int (*FD_MODELINIT)(PARARRAY pars, 
			    FILE * stream, 
			    IMODEL & model);

/** corresponding destructor */
typedef void (*FD_MODELDEST)(void ** specs);


/** Computes time step for internal simulation time grid. This must be
    possible based on information contained in pars, eg. max wave
    velocity, scheme type and order, etc., and primal spatial grid g,
    eg. space steps.
    
    @param[in] pars - parameter array, assumed initialized.
    
    @param[in] stream - verbose output stream
    
    @param[in] g - primal simulation \ref grid , initialized via I/O
    on first IOKEY
    
    @param[out] dt - time step, will have been initialized on call by prior
    call to FD_TIMEGRID
    
    @return - 0 on success, else error code.
    
*/
typedef int (*FD_TIMEGRID)(PARARRAY * pars, 
			   FILE * stream, 
			   grid const & g, 
			   ireal & dt);

  /** creates FD stencils. A \ref stencil describes the dependencies
      between arrays participating in a finite difference scheme, in
      detail. IWAVE uses this information to create ghost cells for
      boundary condition implementation, and for data exchange between
      processors in a domain decomposition.

      See documentation for \ref stencil for a detailed description of stencil construction.

      Preconditions - output stream open, model spatial dimension available,
      grid type and stencil dependency arrays initialized. Follows calls to 
      set_grid_type and build_sten_dep.
      
      Postconditions - \ref stencil object initialized

      @param[in] specs - \ref IMODEL .specs struct, containing
      model-dependent info 

      @param[in] stream - verbose output stream

      @param[in] ndim - model spatial dimension 

      @param[in] gtype - grid type array - assigned from FIELDS - for
      each dim, 0 for primal, 1 for dual (staggered)

      @param[out] sten - stencil object

      @return - 0 for success, else error code

      called in fd_modelcrea
  */
  typedef int (*FD_STENCIL)(void * specs,
			    FILE * stream,
			    int ndim, 
			    IPNT gtype[RDOM_MAX_NARR], 
			    STENCIL * sten);

/** sanity check for coefficient fields - called after
   these are initialized, at beginning of time loop. Should
   throw RVLException for bound transgression or other sin. 
   Parameters for tests stored in model->specs. These 
   tests a priori refer only to the reference RDOM.

   @param dom    - reference RDOM, containing reference simulation fields
   
   @param specs  - \ref IMODEL .specs object

   @param stream - verbose output stream

   called in IWaveSim::run
*/
typedef void (*FD_CHECK)(RDOM * dom,
			 void * specs,
			 FILE * stream);

/*
  IWAVE field spec struct. Model field structure 
  by static array of these. By convention the reference
  grid corresponds to the first of these fields.
*/

typedef struct s_field {
  std::string field;
  int dynamic;
  int substep;
  IPNT gtype;
} FIELD;

/* modeling task tuple - pertains only to forward 
   model. Couples keywords to rarray indices. Task can define 
   input or output, and can be active (coordinates of point in 
   operator domain) or passive (auxiliary parameters needed in
   simulation 

   i/o harness for simulation spec'd via a static array of these.
   by convention the first element should give the input connection
   for the first model field, and provides the keyword for looking
   up the external data source for the primary model grid.
*/
typedef struct s_iokeys {
  std::string keyword;
  int rarrindex;
  bool input;
  bool active;
} IOKEY;

/** basic model definition class
*/
class IWaveInfo {
public:
  /** critical data members; must be defined in 
      glabal namespace somewhere
  */
  static std::string iwave_model;
  static FIELD iwave_fields[]; 
  static IOKEY iwave_iokeys[];
  static FD_MODELINIT minit;
  static FD_MODELDEST mdest;
  static FD_TIMESTEP timestep;
  static FD_TIMEGRID timegrid;
  static FD_STENCIL createstencil;
  static FD_CHECK check;
  // these are redundant relic of previous design, but don't want to chase
  // them down right now - 05.01.14
  std::string get_iwave_model() const { return this->iwave_model; }
  IOKEY * get_iwave_iokeys() const { return &(this->iwave_iokeys[0]); }
  FIELD * get_iwave_fields() const { return &(this->iwave_fields[0]); }
  FD_MODELINIT get_minit() const { return this->minit; }
  FD_MODELDEST get_mdest() const { return this->mdest; }
  FD_TIMEGRID get_timegrid() const { return this->timegrid; }
  FD_TIMESTEP get_timestep() const { return this->timestep; }
  FD_STENCIL get_stencil() const { return this->createstencil;}
  FD_CHECK get_check() const { return this->check; }
  // defined
  int get_num_fields() const;
  int get_num_iokeys() const;
  ostream & write_iwave_fields(ostream & str) const;
  ostream & write_iwave_iokeys(ostream & str) const;
};

namespace TSOpt {

  using RVL::RVLException;
  using RVL::Writeable;

  /* task component: defines relation between IWAVEs in
     IWaveTreeState, keywords indexing i/o info (filenames), input
     flags, and indices of rarrays in the simulation RDOM.
   */
  typedef struct s_task_reln {
    int iwaveindex;   
    std::string keyword;
    int rarrindex;
    bool input;
  } TASK_RELN;

  /* generates list of task relations defining i/o patterns for
     derivatives of a model. IO represented by a IODefn
     object. Algorithm: for order i, loops over 2^i copies of IO.
     In all cases, passive relations unaltered from basic model - all
     belong in order 0 branch. Forward mode: active inputs in comp 0
     (reference model), 1, 2, 4,...,2^{i-1}, with same access pattern
     as basic model. Keyword in comp 2^k altered by addition of suffix
     _d{k}. Active outputs in comp 2^{i}-1. Adjoint mode: inputs same
     except for comp 2^{i-1}, which becomes output with keyword suffix
     _b{i-1}. All forward ouputs in comp 2^{i-1} become inputs.
  */
    
  void IOTask(std::vector<TASK_RELN *> & tr, int order, bool fwd, IWaveInfo const & ic);

  void IOTaskWriter(std::vector<TASK_RELN *> const & tr, ostream & str);

}

#endif
