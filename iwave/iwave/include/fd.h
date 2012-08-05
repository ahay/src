/*
  Author: XW 01-15-2010
  wangxin.tom@gmail.com
  fd.h 
  explicit finite difference framework
*/
/*===========================================================================*/
/**
 * \file fd.h
 * Explicit finite difference modeling framework
 *
 * Contains the \ref FD_MODEL and \ref FD_TS_PARS structs
 */

#ifndef __XWAVE_FD_H_
#define __XWAVE_FD_H_
/*---------------------------------------------------------------------------*/

#include "rdomain.h"
#include "exchange.h"
#include "model.h"
#include "gridio.h"

//#define VERBOSE
/** computational domain includes the boundary points or not */
#define INCLUDE_BOUNDARY_PNTS 0
/** value dependency */
#define DEP_F      4
/** z-derivative dependency */
#define DEP_DFDZ   1
/** x-derivative dependency */
#define DEP_DFDX   2
/** y-derivative dependency */
#define DEP_DFDY   3
  
#define DUAL_GRID   1
#define PRIMAL_GRID 0

/*---------------------*/

/*----------------------------------------------------------------------------*/
/*
  
*/
/**
 Time step function. 
 
 @param[out] dom - RDOM storing dynamic and static fields
 @param[in] iv - internal (sub)step index
 @param[in] fdpars - opaque FD_MODEL parameter object 
 @return 0 on success, else error code

 called in iwave_run. 

 */
typedef int (*TIMESTEP_FUN)(RDOM * dom, int iv, void* fdpars);

/*----------------------------------------------------------------------------*/

/** Model - function members.

 Generic structure for explicit finite difference modeling for 1st order wave equations.
 Provides unassigned variables and function pointers with defined interfaces.
 
 Functions as an abstract base class: construct a concrete model type by implementing all
 function interfaces and assigning function pointers.
 
 Assumptions on computational grid:
 <ul>
 <li>Physical domain's top-left corner is a gridpoint, with indices [0, 0, 0],</li>
 <li>discrete variables can be defined on either the primal grid or one of its dual grids:
 <ul>
 <li>the primal grid points have integer indices,</li>
 <li>the dual grid points have half-integer indices in some dimensions and integer 
 indices in others.</li>
 <li>physical domain boundaries lie on the primal grid</li>
 </ul>

This header file includes a list of functions to be defined, in
cut-paste format. A concrete model type will define a function with
each of these signatures, then assign the corresponding function
pointer from the list to its concrete defined function. These
functions will be referenced only via their assigned pointers in the
FD_MODEL struct, so no header file declarations are necessary.

One of the function interfaces, fd_model_init, acts as a
constructor. The concrete function assigned to this pointer should
assign all function pointers listed in FD_MODEL, including
fd_model_init itself, to defined functions.

*/
/*
 list of functions to be defined, in cut-paste format:
 
 int (*isarr)(int);
 int (*isdyn)(int);
 const char* (*ind2str)(int);  
 int (*set_grid_type)(FILE *, int, IPNT[RDOM_MAX_NARR]);
 int (*build_sten_dep)(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
 int (*create_sten)(FILE *, int, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *, void *);
 int (*readgrid)(PARARRAY *, FILE *, IMODEL *);
 int (*assign_action_array)(FILE *, IMODEL *);
 int (*readmedia)(PARARRAY *, FILE *, IMODEL *, int);
 int (*readtimegrid)(PARARRAY *, FILE *, IMODEL *);
 int (*readschemeinfo)(PARARRAY *, FILE *, struct FD_MODEL *, IMODEL *);
 int (*alter_dom)(int, IPNT, IPNT);
 int (*tsf)(RDOM *,int,FD_MODEL *);
 int (*fd_model_init)(PARARRAY *,FILE *,FD_MODEL *,IMODEL *);
 int (*fd_model_dest)(FD_MODEL *);
*/

/*
 <li>relation between variables is determined by stencil dependence matrix in the way that
 *   if variable ir is dependent of variable ip, sten_dep_mat[ir][ip] is 
 *   assigned a positive number which indicates the dependency type, e.g.,
 *   \ref DEP_F, \ref DEP_DFDZ, \ref DEP_DFDX, \ref DEP_DFDY, etc.
 * - Physical + artificial domain's:
 *    - top-left corner: [-nxl, -nyl, -nzl], if INCLUDE_BOUNDARY_PNTS = 1, 
 *      otherwise [-nxl+1, -nyl+1, -nzl+1]
 *    - bottom-right corner: [nx+nxr-1, ny+nyr-1, nz+nzr-1], if INCLUDE_BOUNDARY_PNTS = 1,
 *      otherwise [nx+nxr-2,ny+nyr-2,nz+nzr-2]
 *    - when INCLUDE_BOUNDARY_PNTS = 0, dynamic variables are not recomputed, but updated using 
 *    boundary conditions
 */
typedef struct FD_MODEL {

  /********* DATA MEMBERS *************/
  /**internal parameters of time step functions - model-dep */
  void * fdpars;

  /********** FUNCTION MEMBERS ************/

  /** returns true (1) if argument is index of array participating in
      simulation, else false (0). Array indices run from 0 to
      RDOM_MAX_NARR - 1 (defined in
      iwave/base/include/utils.h). Typically only a subset of this
      maximal set of indices corresponds to arrays used in a
      simulation. For example, a 2D constant density acoustics
      simulation might use three arrays - pressure at current time
      level (0), pressure at previous time level (1), and sound
      velocity (or square slowness or some other equivalent) (2). Thus
      in this case isarr(i)=1 if i=0, 1, or 2, =0 otherwise. In some
      apps with data-dependent spatial dimension, the number of active
      arrays and their indices may depend on dimension, and possibly
      on other attributes of the simulation available on
      initialization. If that is the case, then fd_model_init should
      determine and store the spatial dimension in a form accessible
      to this function and several others, eg. as static
      data. Therefore calling fd_model_init is made a precondition.
      
      Preconditions: number of active arrays is less than or equal to
      RDOM_MAX_NARR, and fd_model_init has been called.

      Postconditions: none.

      @param[in] i - array index
      @return 1 if i is index of active array, else 0

      Passed as function pointer argument to fd_setcompdom
      (iwave/src/model/fd.c), where it is used in determining which
      array bounds to compute.
  */
  int (*isarr)(int i);

  /** returns string (name or mnemonic) of field indexed by arg in
      rdomain. May depend on spatial dimension.

      Preconditions: fd_model_init has been called 

      Postconditions: none

      @param[in] i - array index
      @return C string (null-terminated char array) if i is index of an active array, else NULL.

      Called in fd_modelcrea (iwave/src/model/fd.c) and sten_out
      (iwave/src/model/stencil.c) to format informative output.
  */
  const char* (*ind2str)(int i);  

  /** returns number of substeps in time step. Presumed to be
      characteristic of timestepping scheme and not dependent on any
      runtime conditions.

      Preconditions: none

      Postconditions: none

      @return number of substeps

      used in iwave_construct to initialize TIMESTEPINFO component of IMODEL
  */
  int (*numsubsteps)();

  /** returns true if array ia is updated in substep iv.

      Preconditions: fd_model_init has been called

      Postconditions: none

      @param ia index of array @param iv substep index 

      @return 1 if array ia is updated in substep iv, 0 otherwise
      (eg. if ia or iv are out of range)

      Called in iwave_run to determine communication patterns. Also isdyn
      is implemented in terms of this function
  */
  int (*update)(int ia, int iv);
      
  /** reads spatial grid info and records it in data members of \ref
      IMODEL struct. IWAVE presumes that all grids involved in the
      simulation are derived from a common primal grid, and this
      function initializes it, normally from a file named in the
      pararray passed as first arg. The function does two jobs:

      (1) Initializes data member IMODEL.g with
      data of primal or reference grid. This object is a <a
      href="../../../grid/doc/html/index.html">grid</a> struct. If
      grid data is stored in an RSF grid file, then the i/o functions
      provided in the <a
      href="../../../grid/doc/html/index.html">grid</a> package can be
      used to do the job.  

      (2) Also initializes padding parameters IMODEL.nls and
      IMODEL.nrs. The primal grid may define the physical domain of
      one of the static fields, without providing extra volume for
      expression of boundary conditions. The padding parameters are
      integer arrays (IPNTs) which store the number of extra grid
      cells along each axis, negative (nls) and positive (nrs)
      directions. These are defaulted to arrays of zeros, so if no
      padding is needed this task can be ignored.

      The grid and padding data may be specified in any convenient
      way. The pararray passed as the first argument, for example, may
      contain a filename for the grid and individual key=value pairs
      for the padding parameters.

      Preconditions: The pararray passed as first argument is
      initialized, and the verbose output stream (second argument) is
      open.

      Postconditions: IMODEL.g, IMODEL.nls, and IMODEL.nrs are initialized.

      @param[in] par - associative array struct containing simulation
      parameters as key=value pairs 

      @param[in] stream - verbose output stream
      
      @param[in] mdl - \ref IMODEL struct, declared at driver level and
      passed as part of IWAVE struct through iwave_construct.

      @return 0 on successful completion, else nonzero error code

      Called in iwave_construct (iwave/src/state/iwave.c). */
  int (*readgrid)(PARARRAY * par, FILE * stream, IMODEL * mdl);

  /** sets grid types for sim fields. Encoded in IPNT (int
      array). Grid may be primal (cell-centered) or dual
      (face-centered) in each axis. Thus a 3D primal field would have
      grid type array [0, 0. 0], whereas a field dual on axis 1 and
      primal on axes 2 and 3 would have grid type array [1, 0, 0]. For
      example, in elastic staggered grid schemes, diagonal stresses
      are primal, velocity on axis i is dual on axis i, primal on
      other axes, and off-diagonal stresses are dual on two axes.

      Preconditions: output stream (1st arg) is open, model dimension
      (2nd arg) is set, grid type array statically allocated in fd_modelcrea.

      Postconditions: grid type array initialized

      @param[in] stream - output stream 
      @param[in] ndim - model dimension (1, 2, or 3)
      @param[out] gtype - array of grid type arrays 

      called in fd_modelcrea
  */
  int (*set_grid_type)(FILE * stream, int ndim, IPNT gtype[RDOM_MAX_NARR]);

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
  int (*build_sten_dep)(FILE * stream, int ndim, int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]);

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
  int (*create_sten)(struct FD_MODEL * fdm,
		     FILE * stream,
		     int ndim, 
		     IPNT gtype[RDOM_MAX_NARR], 
		     int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		     STENCIL * sten);

  /** model-specific alterations of standard domain geom bound to
      grid. The logic in fd_setcompdom can only define arrays based on
      primal and dual grids. Any modifications, such as lowered
      dimensions, must be applied after all computational arrays are
      declared in fd_setcompdom, and exchange buffers set up in
      ex_compute, but before memory is allocated in fd_modelcrea. This
      function provides an interface for any such array
      modifications. The first argument is the index of the array to
      be (possibly) modified, the second and third arguments are the
      start and end index arrays for this array.

      For example, PML coefficient arrays may be effectively
      one-dimensional. This routine may be used to make them actually
      one-dimensional, thus reducing storage substantially. For the
      PML coefficient regulating decay along the second axis, one
      might write

      gs[0]=gs[1]; ge[0]=ge[1];
      for (i=1;i<RARR_MAX_NDIM;i++) { gs[i]=0; ge[i]=0; },

      thus transferring the start and end indices of the second axis
      in a multidimensional array to the first axis and trivializing
      the other axes - after which the array is actually
      one-dimensional.

      Preconditions: all nontrivial data-dependent info is already
      encoded in gs and ge - effectively, this function can only
      reorder the information already encoded or zero some of it. Thus
      fd_setcompdom and ex_compute must already have been called.

      Postconditions: start, end indices for array ia are correct,
      memory may be allocated.

      @param[in] ia - index of array to be adjusted
      @param[out] gs - start indices of array
      @param[out] ge - end indices of array
      @return - 0 on success, else error code

      Note that if no arrays need alteration, this function may be
      implemented as a no-op. If there were a base class, it would be
      a virtual function with no-op implementaiton which could be
      overrridden.

      called in fd_modelcrea
 */
  int (*alter_dom)(int ia, IPNT gs, IPNT ge);

  /** initializes static arrays, typically coefficients in the system
      of PDEs being solved. Data to intialize these arrays may be read
      from files, or computed from auxiliary information stored in the
      \ref PARARRAY argument pars. The static arrays are \ref RARR s stored
      in the \ref IMODEL :: \ref RDOM ld_a (allocated domain). This
      interface is very flexible and can accommodate many methods for
      initializing static data.

      Current IWAVE applications access detailed array data via the
      RSF file format. The I/O functions provided in the <a
      href="../../../grid/doc/html/index.html">grid</a> package assign
      values to array entries when these correspond to gridpoints
      represented in the file data, and extend the data systematically
      to array entries outside of the file-defined grid. Thus the
      array need not be the same size or shape as the grid - it is
      however assumed that the sample rate implicitly represented in
      the array is that of the grid.

      The implementation of this function is responsible for details
      such as interpolation to half-integer gridpoints for dual grid
      arrays (if necessary).

      TODO: interpolation could be flagged by passing the gtype
      array to the gridio functions, and coded generically, which
      would remove a lot of complexity from implementations of
      readmedia.

      Preconditions: PARARRAY argument initialized, stream open. Data
      structures in the \ref IMODEL arg must be allocated, so
      readgrid and several other functions must have been called
      first. In particular the \ref grid object mdl.g must be
      initialized, and the RDOM mdl.ld_a as well.

      Postconditions: static arrays in mdl.ld_a initialized.

      @param[in] pars - parameter array, assumed initialized. Source of
      filenames, for example - if "velocity" names one of the static
      arrays, then the pair velocity=filename should signal that the
      "velocity" array is to be initialized from the file with name
      filename. Could also include constant values, or polynomial
      coefficients, or any other convenient encoding of static array
      values.

      @param[in] stream - verbose output stream

      @param[out] mdl - \ref IMODEL whose static arrays are being initialized

      @param[in] panelnum - for multisimulations in which different
      models are used in different simulations, this parameter should
      indicate dhte simulation index. For ordinary simulations
      depending on one model, default = 0.

      @return - 0 on success, else error code.

      called in fd_mread 
  */
  int (*readmedia)(PARARRAY * pars, FILE * stream, IMODEL * mdl, int panelnum);

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

      called in fd_mread
  */
  int (*readtimegrid)(PARARRAY * pars, FILE * stream, IMODEL * mdl);

  /** fills in scheme info dependent on spatial grid, time step and
      coefficient arrays. This information (for example, Courant
      ratios) may be stored in the opaque data member FDPARS of
      FD_MODEL.

      Preconditions: all geometry, static array, and time step info
      has been initialized - must be called after readtimegrid.

      Postconditions: any scheme dependent info not already
      initialized is set - final step in scheme initializiation.

      @param[in] pars - parameter array, assumed initialized.

      @param[in] stream - verbose output stream

      @param[out] mdl - \ref IMODEL whose static arrays are being initialized
      
      @return - 0 on success, else error code.

      called in fd_mread
  */
  int (*readschemeinfo)(PARARRAY * pars, FILE * stream, IMODEL * mdl);

  /** timestep function - called in iwave_run */
  TIMESTEP_FUN tsf;

  /** parameter copy function - since opaque FDPARS parameter struct is unrestricted,
      a function to define copying is necessary */
  void (*parcopy)(void *, const void *);

  /** main constructor - should do at least these three things: - 

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
      
      @return - 0 on success, else error code.
  */
  int (*fd_model_init)(PARARRAY * pars, 
		       FILE * stream, 
		       IMODEL * mdl);

  /** destructor - also virtual, since fd models may allocate memory
      (notably for fdpars) 
  */
  int (*fd_model_dest)(struct IMODEL *);
    
} FD_MODEL;

/** default constructor (implemented)*/
void fd_model_setnull(FD_MODEL * fd);

/** (implemented) returns true if arg is index in rdomain of a dynamic field
    (i.e. a field updated in the simulation), else false. For
    example, in the constant density acoustic example, the current
    and previous time levels of pressure are dynamic, whereas the
    sound velocity is not. So isdyn(i)=1 if i=0 or 1, 0 else. May
    depend on spatial dimension, so same preconditions as isarr.
    
    Preconditions: fd_model_init has been called
    
    Postconditions: none

    @param[in] fd - FD_MODEL defining scheme
    @param[in] i - array index
    @return 1 if i is index of dynamic field, else 0
    
    Called in fd_modelcrea (only dynamic arrays need ghost cell
    augmentation) and fprint_weqn (iwave/src/model/fd.c)
*/
int isdyn(FD_MODEL * fd, int i);

/** number of dynamic arrays defined in scheme

    @param[in] fd - FD_MODEL defining scheme
    @return number of dynamic arrays

Called in setrecvexchange (iwaveinfo.c)
*/
int narr(FD_MODEL * fd);

/*--------------*/
#endif
