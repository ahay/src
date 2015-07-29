#ifndef __IWAVE_M_TUPLE
#define __IWAVE_M_TUPLE

#include "std_cpp_includes.hh"
#include "utils.h"
#include "except.hh"
#include "write.hh"
#include "model.h"

using RVL::RVLException;

size_t pow2(int);

/* max number of data tuples permitted */
#define IWAVEMAXDATA 256

// forward declaration
class IWaveInfo;

/**
 Time step function. 
 
 @param[out] dom - std::vector<RDOM *> storing dynamic and static fields
 @param[in] iv - internal (sub)step index
 @param[in] fdpars - opaque FD_MODEL parameter object 
 @return 0 on success, else error code

 called in iwave_run. 

 */
typedef void (*FD_TIMESTEP)(std::vector<RDOM *> dom, bool fwd, int iv, void* fdpars);

/** FD model constructor - should do at least these three things: - 
    
    <ul> 
    
    <li>assign FD_MODEL ("this") to IMODEL.specs (would do this in a
    base class if there were such a thing)</li>
    
    <li> assign all function pointers in FD_MODEL, including this
    one </li>
    
    <li> create appropriate data structure to serve as fdpars
    parameter struct for FD_MODEL specialization, initialize data
    members of FD_MODEL.fdpars using data from pars, the IMODEL data
    members, or any other source, as required.</li> 
    
    </ul>
    
    @param[in] pars - parameter array, assumed initialized.
    
    @param[in] stream - verbose output stream
    
    @param[out] mdl - \ref IMODEL whose static arrays are being initialized
    
    @param[in] gridkey - key string for source of grid model info
    
    @return - 0 on success, else error code.
*/

typedef int (*FD_MODELINIT)(PARARRAY * pars, 
			    FILE * stream, 
			    grid const & g,
			    ireal dt,
			    void ** specs);

/** corresponding destructor */
typedef void (*FD_MODELDEST)(void ** specs);


  /** initializes \ref TIMEINDEX data member of \ref IMODEL mdl - this
      chiefly means compute time step. Might be read from pars, or
      computed from velocity or other fields initialized by readmedia,
      hence called after readmedia in fd_mread.

      Preconditions: output stream open, \ref PARARRAY pars
      initialized, \ref IMODEL mdl initialized.

      Postconditions: mdl.tsind.dt and other data members set
      appropriately.

      @param[in] pars - parameter array, assumed initialized.

      @param[in] stream - verbose output stream

      @param[out] mdl - \ref IMODEL whose static arrays are being initialized
      
      @return - 0 on success, else error code.

  */
typedef int (*FD_TIMEGRID)(PARARRAY * pars, 
			   FILE * stream, 
			   grid const & g, 
			   ireal & dt);

  /** assigns stencil dependency array. This array is auxiliary
      information used in creating a stencil struct. It describes the
      relations between fields, at the level of the PDEs. Since IWAVE
      applies to linear first-order systems of pdes, each field can
      depend on other fields either through their values or through
      the values of their partial derivatives. These dependency types
      are encoded as follows:

      <ul>
      <li> dependent on value -> DEP_F</li>
      <li> dependent on 1st derivative wrt axis-z -> DEP_DFDZ</li>
      <li> dependent on 1st derivative wrt axis-x -> DEP_DFDX</li>
      <li> dependent on 1st derivative wrt axis-y -> DEP_DFDY</li>
      </ul>

      The [i][j] element of the dependency matrix (output) contains
      the type code for the dependence of field i on field j. The type
      codes DEP_* are defined in this file (fd.h).

      Preconditions: output stream open, model dimension set,
      dependency matrix statically allocated in fd_modelcrea.

      Postconditions: dependency matrix initialized

      @param[in] stream - output stream 
      @param[in] ndim - model dimension (1, 2, or 3)
      @param[out] stendep - stencil dependency matrix
      @return - 0 for success, else error code.

      called in fd_modelcrea
  */
//  int (*build_sten_dep)(FILE * stream, int ndim, int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]);

  /** creates FD stencils. A \ref stencil describes the dependencies
      between arrays participating in a finite difference scheme, in
      detail. IWAVE uses this information to create ghost cells for
      boundary condition implementation, and for data exchange between
      processors in a domain decomposition.

      Preconditions - output stream open, model spatial dimension available,
      grid type and stencil dependency arrays initialized. Follows calls to 
      set_grid_type and build_sten_dep.
      
      Postconditions - \ref stencil object initialized

      @param[in] stream - verbose output stream
      @param[in] ndim - model spatial dimension
      @param[in] gtype - grid type array
      @param[in] stendep - stencil dependency matrix 
      @param[out] sten - stencil object
      @param[in] fdpars - opaque FD parameter object
      @return - 0 for success, else error code

      called in fd_modelcrea
  */
  typedef int (*FD_STENCIL)(void * specs,
			    FILE * stream,
			    IWaveInfo const & ic,
			    int ndim, 
			    IPNT gtype[RDOM_MAX_NARR], 
			    //		     int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR], 
			    STENCIL * sten);

/* sanity check for coefficient fields - called after
   these are initialized, at beginning of time loop. Throws 
   suitable exception for bound transgression or other sin. 
   Parameters for tests stored in model->specs. These 
   tests a priori refer only to the reference RDOM.
*/
typedef void (*FD_CHECK)(RDOM * dom,
			 void * specs,
			 FILE * sream);

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
