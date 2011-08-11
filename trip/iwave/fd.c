/* Explicit finite difference modeling framework */

/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/**
 * XW 01-18-2010
 * wangxin.tom@gmail.com
 * fd.c
 */

#include "fd.h"

#include "model.h"
#include "stencil.h"
/*^*/

#include "exchange.h"

#ifndef _sf_fd_h

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
/*^*/

/**
   Time step function.
   *
   @param[out] pdom - RDOM in which to store output (updated) dynamic fields (next time level)
   @param[in] rdom - RDOM from which to take input dynamic fields (current time level)
   @param[in] cdom - RDOM from which to take input static fields (coefficient arrays)
   @param[in] ia - index of array to be updated
   @param[in] fdpars - opaque FD_MODEL parameter object 
   @return 0 on success, else error code
   *
   called in iwave_run
*/
typedef int (*TIMESTEP_FUN)(RDOM * pdom, RDOM * rdom,RDOM * cdom,int ia, void* fdpars);
/*^*/

/**
   Post time step function
   *
   Time step functions update the computational \ref RDOM. If not
   identical to the allocated \ref RDOM, then another function is
   required to update the latter, typically by implementing boundary
   conditions on the ghost cells which form the set-theoretic
   difference. Both allocated and computational domains, plus much
   other information about the model, may figure into this
   computation, so the argument is the entire \ref IMODEL.
   *
   @param ia - index of array to be updated
   @param mdl - simulation \ref IMODEL
   @return - 0 on success, else error code
   *
   called in iwave_run
*/
typedef int (*POSTTIMESTEP_FUN)(int ia, IMODEL * mdl);
/*^*/

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
	*
	Preconditions: number of active arrays is less than or equal to
	RDOM_MAX_NARR, and fd_model_init has been called.
	*
	Postconditions: none.
	*
	@param[in] i - array index
	@return 1 if i is index of active array, else 0
	*
	Passed as function pointer argument to fd_setcompdom
	(iwave/src/model/fd.c), where it is used in determining which
	array bounds to compute.
    */
    int (*isarr)(int i);
    /** returns true if arg is index in rdomain of a dynamic field
	(i.e. a field updated in the simulation), else false. For
	example, in the constant density acoustic example, the current
	and previous time levels of pressure are dynamic, whereas the
	sound velocity is not. So isdyn(i)=1 if i=0 or 1, 0 else. May
	depend on spatial dimension, so same preconditions as isarr.
	*
	Preconditions: fd_model_init has been called
	*
	Postconditions: none
	*
	@param[in] i - array index
	@return 1 if i is index of dynamic field, else 0
	*
	Called in fd_modelcrea (only dynamic arrays need ghost cell
	augmentation) and fprint_weqn (iwave/src/model/fd.c)
    */
    int (*isdyn)(int i);
    /** returns string (name or mnemonic) of field indexed by arg in
	rdomain. May depend on spatial dimension.
	*
	Preconditions: fd_model_init has been called 
	*
	Postconditions: none
	*
	@param[in] i - array index
	@return C string (null-terminated char array) if i is index of an active array, else NULL.
	*
	Called in fd_modelcrea (iwave/src/model/fd.c) and sten_out
	(iwave/src/model/stencil.c) to format informative output.
    */
    const char* (*ind2str)(int i);  
    /** reads spatial grid info and records it in data members of \ref
	IMODEL struct. IWAVE presumes that all grids involved in the
	simulation are derived from a common primal grid, and this
	function initializes it, normally from a file named in the
	pararray passed as first arg. The function does two jobs:
	*
	(1) Initializes data member IMODEL.g with
	data of primal or reference grid. This object is a <a
	href="../../../grid/doc/html/index.html">grid</a> struct. If
	grid data is stored in an RSF grid file, then the i/o functions
	provided in the <a
	href="../../../grid/doc/html/index.html">grid</a> package can be
	used to do the job.  
	*
	(2) Also initializes padding parameters IMODEL.nls and
	IMODEL.nrs. The primal grid may define the physical domain of
	one of the static fields, without providing extra volume for
	expression of boundary conditions. The padding parameters are
	integer arrays (IPNTs) which store the number of extra grid
	cells along each axis, negative (nls) and positive (nrs)
	directions. These are defaulted to arrays of zeros, so if no
	padding is needed this task can be ignored.
	*
	The grid and padding data may be specified in any convenient
	way. The pararray passed as the first argument, for example, may
	contain a filename for the grid and individual key=value pairs
	for the padding parameters.
	*
	Preconditions: The pararray passed as first argument is
	initialized, and the verbose output stream (second argument) is
	open.
	*
	Postconditions: IMODEL.g, IMODEL.nls, and IMODEL.nrs are initialized.
	*
	@param[in] par - associative array struct containing simulation
	parameters as key=value pairs 
	*
	@param[in] stream - verbose output stream
	*
	@param[in] mdl - \ref IMODEL struct, declared at driver level and
	passed as part of IWAVE struct through iwave_construct.
	*
	@return 0 on successful completion, else nonzero error code
	*
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
    *
    Preconditions: output stream (1st arg) is open, model dimension
    (2nd arg) is set, grid type array statically allocated in fd_modelcrea.
    *
    Postconditions: grid type array initialized
    *
    @param[in] stream - output stream 
    @param[in] ndim - model dimension (1, 2, or 3)
    @param[out] gtype - array of grid type arrays 
    *
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
    *
    <ul>
    <li> dependent on value -> DEP_F</li>
    <li> dependent on 1st derivative wrt axis-z -> DEP_DFDZ</li>
    <li> dependent on 1st derivative wrt axis-x -> DEP_DFDX</li>
    <li> dependent on 1st derivative wrt axis-y -> DEP_DFDY</li>
    </ul>
    *
    The [i][j] element of the dependency matrix (output) contains
    the type code for the dependence of field i on field j. The type
    codes DEP_* are defined in this file (fd.h).
    *
    Preconditions: output stream open, model dimension set,
    dependency matrix statically allocated in fd_modelcrea.
    *
    Postconditions: dependency matrix initialized
    *
    @param[in] stream - output stream 
    @param[in] ndim - model dimension (1, 2, or 3)
    @param[out] stendep - stencil dependency matrix
    @return - 0 for success, else error code.
    *
    called in fd_modelcrea
*/
    int (*build_sten_dep)(FILE * stream, int ndim, int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]);
/** creates FD stencils. A \ref stencil describes the dependencies
    between arrays participating in a finite difference scheme, in
    detail. IWAVE uses this information to create ghost cells for
    boundary condition implementation, and for data exchange between
    processors in a domain decomposition.
    *
    Preconditions - output stream open, model spatial dimension available,
    grid type and stencil dependency arrays initialized. Follows calls to 
    set_grid_type and build_sten_dep.
    * 
    Postconditions - \ref stencil object initialized
    *
    @param[in] stream - verbose output stream
    @param[in] ndim - model spatial dimension
    @param[in] gtype - grid type array
    @param[in] stendep - stencil dependency matrix 
    @param[out] sten - stencil object
    @param[in] fdpars - opaque FD parameter object
    @return - 0 for success, else error code
    *
    called in fd_modelcrea
*/
    int (*create_sten)(FILE * stream,
		       int ndim, 
		       IPNT gtype[RDOM_MAX_NARR], 
		       int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		       STENCIL * sten, void * fdpars);
/** initializes the \ref TIMESTEPINFO data member of \ref
    IMODEL. The data members of \ref TIMESTEPINFO are:
    *
    <ul>
    <li>int narr - number of dynamic arrays</li>
    *
    <li>int arrs[RDOM_MAX_NARR] - list of dynamic array indices -
    only first narr entries are significant</li>
    *
    <li>int npair - number of action pairs, corresponding to
    distinct actions during timestep loop</li>
    *
    <li>TIMESTEPINFOPAIR pairs[MAX_ACTIONS] array of action pairs,
    defining simulation timestep </li>
    *
    </ul>
    *
    [NOTE: narr and arrs are actually redundant, and could be
    replaced by computations involving the isdyn function. Requires
    refactoring iwavefun.c and iwaveinfo.c.]
    *
    A \ref TIMESTEPINFOPAIR consists of two int data members, arr
    and action. The arr member is the index of a dynamic array. The
    action member is an action code, selected from those defined in
    \ref model_actions.h . The present version of IWAVE defines only
    two actions, ACTION_COMPUTE and ACTION_EXCHANGE.
    *
    Example: if a simulation involves two dynamic arrays, a and b,
    with indices 0 and 1, and each is updated in each domain then
    exchanged with neighboring domains, then npair=4 and
    *
    <ol>
    <li>update a - pairs[0].arr=0, pairs[0].action=ACTION_COMPUTE</li>
    <li>exchange a - pairs[1].arr=0, pairs[1].action=ACTION_EXCHANGE</li>
    <li>update b - pairs[2].arr=1, pairs[2].action=ACTION_COMPUTE</li>
    <li>exchange b - pairs[3].arr=1, pairs[3].action=ACTION_EXCHANGE</li>
    </ol>
    *
    Preconditions: output stream open
    *
    Postconditions: TIMESTEPINFO object initialized
    *
    @param stream[in] - verbose output stream
    * 
    @param im[out] - \ref IMODEL struct containing \ref TIMESTEPINFO
    data member to be initialized
    *
    @return - 0 on success, else error code
    *
    called in fd_modelcrea 
*/
    int (*assign_action_array)(FILE * stream, IMODEL * im);
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
    *
    For example, PML coefficient arrays may be effectively
    one-dimensional. This routine may be used to make them actually
    one-dimensional, thus reducing storage substantially. For the
    PML coefficient regulating decay along the second axis, one
    might write
    *
    gs[0]=gs[1]; ge[0]=ge[1];
    for (i=1;i<RARR_MAX_NDIM;i++) { gs[i]=0; ge[i]=0; },
    *
    thus transferring the start and end indices of the second axis
    in a multidimensional array to the first axis and trivializing
    the other axes - after which the array is actually
    one-dimensional.
    *
    Preconditions: all nontrivial data-dependent info is already
    encoded in gs and ge - effectively, this function can only
    reorder the information already encoded or zero some of it. Thus
    fd_setcompdom and ex_compute must already have been called.
    *
    Postconditions: start, end indices for array ia are correct,
    memory may be allocated.
    *
    @param[in] ia - index of array to be adjusted
    @param[out] gs - start indices of array
    @param[out] ge - end indices of array
    @return - 0 on success, else error code
    *
    Note that if no arrays need alteration, this function may be
    implemented as a no-op. If there were a base class, it would be
    a virtual function with no-op implementaiton which could be
    overrridden.
    *
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
    *
    Current IWAVE applications access detailed array data via the
    RSF file format. The I/O functions provided in the <a
    href="../../../grid/doc/html/index.html">grid</a> package assign
    values to array entries when these correspond to gridpoints
    represented in the file data, and extend the data systematically
    to array entries outside of the file-defined grid. Thus the
    array need not be the same size or shape as the grid - it is
    however assumed that the sample rate implicitly represented in
    the array is that of the grid.
    *
    The implementation of this function is responsible for details
    such as interpolation to half-integer gridpoints for dual grid
    arrays (if necessary).
    *
    TODO: interpolation could be flagged by passing the gtype
    array to the gridio functions, and coded generically, which
    would remove a lot of complexity from implementations of
    readmedia.
    *
    Preconditions: PARARRAY argument initialized, stream open. Data
    structures in the \ref IMODEL arg must be allocated, so
    readgrid and several other functions must have been called
    first. In particular the \ref grid object mdl.g must be
    initialized, and the RDOM mdl.ld_a as well.
    *
    Postconditions: static arrays in mdl.ld_a initialized.
    *
    @param[in] pars - parameter array, assumed initialized. Source of
    filenames, for example - if "velocity" names one of the static
    arrays, then the pair velocity=filename should signal that the
    "velocity" array is to be initialized from the file with name
    filename. Could also include constant values, or polynomial
    coefficients, or any other convenient encoding of static array
    values.
    *
    @param[in] stream - verbose output stream
    *
    @param[out] mdl - \ref IMODEL whose static arrays are being initialized
    *
    @param[in] panelnum - for multisimulations in which different
    models are used in different simulations, this parameter should
    indicate dhte simulation index. For ordinary simulations
    depending on one model, default = 0.
    *
    @return - 0 on success, else error code.
    *
    called in fd_mread 
*/
    int (*readmedia)(PARARRAY * pars, FILE * stream, IMODEL * mdl, int panelnum);
/** initializes \ref TIMEINDEX data member of \ref IMODEL mdl - this
    chiefly means compute time step. Might be read from pars, or
    computed from velocity or other fields initialized by readmedia,
    hence called after readmedia in fd_mread.
    *
    Preconditions: output stream open, \ref PARARRAY pars
    initialized, \ref IMODEL mdl initialized.
    *
    Postconditions: mdl.tsind.dt and other data members set
    appropriately.
    *
    @param[in] pars - parameter array, assumed initialized.
    *
    @param[in] stream - verbose output stream
    *
    @param[out] mdl - \ref IMODEL whose static arrays are being initialized
    * 
    @return - 0 on success, else error code.
    *
    called in fd_mread
*/
    int (*readtimegrid)(PARARRAY * pars, FILE * stream, IMODEL * mdl);
/** fills in scheme info dependent on spatial grid, time step and
    coefficient arrays. This information (for example, Courant
    ratios) may be stored in the opaque data member FDPARS of
    FD_MODEL.
    *
    Preconditions: all geometry, static array, and time step info
    has been initialized - must be called after readtimegrid.
    *
    Postconditions: any scheme dependent info not already
    initialized is set - final step in scheme initializiation.
    *
    @param[in] pars - parameter array, assumed initialized.
    *
    @param[in] stream - verbose output stream
    *
    @param[out] mdl - \ref IMODEL whose static arrays are being initialized
    * 
    @return - 0 on success, else error code.
    *
    called in fd_mread
*/
    int (*readschemeinfo)(PARARRAY * pars, FILE * stream, IMODEL * mdl);
/** timestep function - conforms to \ref TIMESTEP_FUN interface,
    called in iwave_run */
    TIMESTEP_FUN tsf;
/** adjoint timestep function - also conforms to \ref TIMESTEP_FUN
    interface. Optional: used only in IWAVE++ inversion extension,
    can be set to NULL in pure modeling applications. */
    TIMESTEP_FUN tsa;
/** post-timestep function - see \ref POSTTIMESTEP_FUN
    typedef. Called in iwave_run */
    POSTTIMESTEP_FUN ptsf;
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
/*^*/

#endif

int fd_setcompdom(FILE * stream, IPNT cdims, IPNT crank, 
                  IMODEL * model, IPNT dgs[], IPNT dge[],
		  //		  int m_size,
		  IPNT gtype[RDOM_MAX_NARR],
		  int (*isarr)(int)) 
/* helper functions */
{
    /**
     * virtual start and end: physical + artificial,  
     * assume physical start at IPNT_0
     */
    IPNT gs_pa, ge_pa;  
    IPNT gn_pa;
    IPNT ls, le; 
    IPNT ns;
    //  int iv, i, idim;
    int i, idim;
    int ndim = (model->g).dim;
  
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	IASN(dgs[i], IPNT_1);
	IASN(dge[i], IPNT_0);
    }
 
    get_n(ns, model->g);
    IASN(gs_pa, IPNT_1);
    IASN(ge_pa, IPNT_0);
    IASN(ls, IPNT_1);
    IASN(le, IPNT_0);
#if INCLUDE_BOUNDARY_PNTS
    /* include left and right boundary pnts */
    /* ls, le: local start, end for primal grid */
    for (idim = 0;idim < ndim;idim ++) {
	gs_pa[idim] = -model->nls[idim];
	ge_pa[idim] = ns[idim] + model->nrs[idim] - 1;
	gn_pa[idim] = ns[idim] + model->nls[idim] + model->nrs[idim];
	ls[idim] = gn_pa[idim] * (long)crank[idim];
	ls[idim] = gs_pa[idim] + ls[idim] / (long)cdims[idim];
	le[idim] = gn_pa[idim] * (long)(crank[idim] + 1);
	le[idim] = gs_pa[idim] + le[idim] / (long)cdims[idim]-1;
	if (le[idim] < ls[idim]) {
	    fprintf(stream, "Error: in fd_setcompdom: le[%d] < ls[%d]\n",
		    idim, idim);
	    return E_INTERNAL;
	}
    } 

    //  for (i = 0;i < m_size;i ++) {
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    iv = getindices(i);
	//    if (!isarr(iv)) continue;
	if (!isarr(i)) continue;
	for (idim = 0;idim < ndim;idim ++) {
	    dgs[i][idim] = ls[idim];
	    dge[i][idim] = le[idim];
	    if ( crank[idim] == cdims[idim]-1 && 
		 gtype[i][idim] == DUAL_GRID )
		dge[i][idim] --;
	}
    }
#else
    /* not include left and right boundary pnts */
    /* ls, le: local start, end for primal grid */
    for (idim = 0;idim < ndim;idim ++) {
	gs_pa[idim] = -model->nls[idim]+1;
	ge_pa[idim] = ns[idim] + model->nrs[idim] - 2;
	gn_pa[idim] = ns[idim] + model->nls[idim] + model->nrs[idim] - 2;
	ls[idim] = gn_pa[idim] * (long)crank[idim];
	ls[idim] = gs_pa[idim] + ls[idim] / (long)cdims[idim];
	le[idim] = gn_pa[idim] * (long)(crank[idim] + 1);
	le[idim] = gs_pa[idim] + le[idim] / (long)cdims[idim]-1;
	if (le[idim] < ls[idim]) {
	    fprintf(stream, "Error: in fd_setcompdom: le[%d] < ls[%d]\n",
		    idim, idim);
	    return E_INTERNAL;
	}
    } 

    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    iv = getindices(i);
	//    fprintf(stream,"fd_setcompdom - array %d index %d\n",i,iv);
	//    if (isarr(iv)) {
	if (isarr(i)) {
	    //      fprintf(stream,"fd_setcompdom - assigning array %d ",iv);
	    //      fprintf(stream,"fd_setcompdom - assigning array %d ",i);
	    for (idim = 0;idim < ndim;idim ++) {
		//	dgs[iv][idim] = ls[idim];
		dgs[i][idim] = ls[idim];
		//	dge[iv][idim] = le[idim];
		dge[i][idim] = le[idim];
		if (crank[idim] == 0 && 
		    //	    gtype[iv][idim] == DUAL_GRID)
		    gtype[i][idim] == DUAL_GRID)
		    //	  dgs[iv][idim] --;
		    dgs[i][idim] --;
		//	fprintf(stream,"dgs[%d]=%d, dge[%d]=%d",idim,dgs[iv][idim],idim,dge[iv][idim]);
		//	fprintf(stream,"dgs[%d]=%d, dge[%d]=%d",idim,dgs[i][idim],idim,dge[i][idim]);
	    }
	    fprintf(stream,"\n");
	}   
    }
#endif

    return 0;
}

//int fprint_weqn(FILE * stream, int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], int (*isdyn)(int), char * (*getnames)(int), int (*getindices)(int), int m_size) {  
int fprint_weqn(FILE * stream, 
		int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		int (*isdyn)(int), 
		char * (*ind2str)(int)) {
    //char * (*getnames)(int), 
    //		int (*getindices)(int),
    //int m_size) {  

    /* try to print the equations */

    //  int i, j, ir, ip;
    int i, j;
    int op;

    fprintf(stream, "\nWAVE EQUATIONS:\n");
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    ir = getindices(i);
	//    if (!isdyn(ir)) continue;
	if (!isdyn(i)) continue;
	fprintf(stream, "d %4s/ dt    =    ", ind2str(i));
	for (j = 0;j < RDOM_MAX_NARR;j ++) {
	    //      ip = getindices(j);
	    //      op = sten_dep_mat[ir][ip];
	    op = sten_dep_mat[i][j];
	    switch (op) {
		case DEP_F:
		    fprintf(stream, "%10s    ", ind2str(j));
		    break;
		case DEP_DFDZ:
		    fprintf(stream, "d %4s/ dz    ", ind2str(j));
		    break;
		case DEP_DFDX:
		    fprintf(stream, "d %4s/ dx    ", ind2str(j));
		    break;
		case DEP_DFDY:
		    fprintf(stream, "d %4s/ dy    ", ind2str(j));
		    break;
		default:
		    break;
	    }
	}
	fprintf(stream,"\n");
    }  
    return 0;
}


void fd_model_setnull(FD_MODEL * fd) 
/*< default constructor >*/
{
    fd->fdpars=NULL;
    /* leave other (statically allocated) data members uninitialised */
    fd->readgrid=NULL;
    fd->readtimegrid=NULL;
    fd->readschemeinfo=NULL;
    fd->set_grid_type=NULL;
    fd->build_sten_dep=NULL;
    fd->create_sten=NULL;
    fd->assign_action_array=NULL;
    fd->readmedia=NULL;
    //  fd->getndim=NULL;
    fd->isdyn=NULL;
    fd->isarr=NULL;
    fd->ind2str=NULL;
    fd->parcopy=NULL;
    fd->fd_model_init=NULL;
    fd->fd_model_dest=NULL;
}

int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model) {
    FD_MODEL *fdm = (FD_MODEL *)model->specs;
    int err;
    int ndim, nnei, inei, idim, iv, i;
    IPNT ns;
    IPNT dgs[RDOM_MAX_NARR], dge[RDOM_MAX_NARR];    /* computational domain */
    IPNT dgsa[RDOM_MAX_NARR], dgea[RDOM_MAX_NARR];  /* allocated domain */
    IPNT gs, ge;
    IPNT dgsr[IWAVE_NNEI], dger[IWAVE_NNEI];                        /* receive domains */
    IPNT dgsrs[IWAVE_NNEI][RDOM_MAX_NARR], dgers[IWAVE_NNEI][RDOM_MAX_NARR];   /* all receive domains */
    int frcvempty[IWAVE_NNEI], rcvne;                       /* empty receive flag */
  
    /* grid type in each dimension: primal grid (=0) or dual grid (=1) */
    IPNT gtype[RDOM_MAX_NARR];         
    /** 
     * stencil dependence matrix:
     * if variable ir is dependent of variable ip, sten_dep_mat[ir][ip] is 
     * assigned a positive number which indicates the dependency type, e.g.,
     * DEP_F, DEP_DFDZ,DEP_DFDX,DEP_DFDY, etc.
     */
    int  sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR];
    /** Stencil - defines size and shape of all FD timestep updates */
    STENCIL sten;
  
    /*--------------------------------------------------------------------------*/
    /*-assign indices for send and receive domains------------------------------*/
    for (iv = 0;iv < IWAVE_NNEI;iv ++) {
	IASN(dgsr[iv], IPNT_1);
	IASN(dger[iv], IPNT_0);
	frcvempty[iv] = 0;
	for (i = 0;i < RDOM_MAX_NARR;i ++) {
	    IASN(dgsrs[iv][i], IPNT_1);
	    IASN(dgers[iv][i], IPNT_0);
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-read grid info-----------------------------------------------------------*/
    /* moved to iwave_construct - 28.01.11 WWS 
       if ( (err=fdm->readgrid(par, stream, model)) ) {
       return err;
       }
    */
    ndim = model->g.dim;

    //  fprintf(stderr,"in fd_modelcrea: ndim=%d\n",ndim);

    /*--------------------------------------------------------------------------*/
    /*-set nnei (num of neighbors in cart grid --*/
    if ( (err=im_setndim(model)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from im_setndim, err=%d\n",err);
	fflush(stream);
	return err;
    }
    nnei = model->nnei;
    /*--------------------------------------------------------------------------*/
    /*-assigning action array---------------------------------------------------*/
    if ( (err=fdm->assign_action_array(stream, model)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from assign_action_array, err=%d\n",err);
	fflush(stream);
	return err;
    }
    /*--------------------------------------------------------------------------*/
    /*-create stencil-----------------------------------------------------------*/
    if ( (err=fdm->set_grid_type(stream,ndim,gtype)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from set_grid_type, err=%d\n",err);
	fflush(stream);
	return err;
    }

    if ( (err=fdm->build_sten_dep(stream,ndim,sten_dep_mat)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from build_sten_dep, err=%d\n",err);
	fflush(stream);
	return err;
    }
  
    //  if ( (err=fdm->create_sten(stream,ndim,fdm->getnarr(),gtype,sten_dep_mat,&sten,fdm->fdpars)) ) {
    if ( (err=fdm->create_sten(stream,ndim,gtype,sten_dep_mat,&sten,fdm->fdpars)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from create_sten, err=%d\n",err);
	fflush(stream);
	return err;
    }

    /* print out stencil if desired */
#ifdef VERBOSE
    sten_out(sten, stream, fdm->ind2str);
    //  fprint_weqn(stream,sten_dep_mat,fdm->isdyn,fdm->getnames,fdm->getindices,fdm->getnarr());
    //  fprint_weqn(stream,sten_dep_mat,fdm->isdyn,fdm->getnames,fdm->getnarr());
    fprint_weqn(stream,sten_dep_mat,fdm->isdyn,fdm->ind2str);
#endif

    /*--------------------------------------------------------------------------*/
    /*-compute local computational grid size------------------------------------*/
    for ( idim = 0; idim < RDOM_MAX_NARR; ++idim ) {
	IASN(dgs[idim], IPNT_1);
	IASN(dge[idim], IPNT_0);
    }
    //  fprintf(stream,"**** in fd_modelcrea: narr=%d\n",fdm->getnarr());
    //  if ( (err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, fdm->getnarr(), gtype,fdm->isarr,fdm->getindices)) ) {
    //  if ( (err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, fdm->getnarr(), gtype,fdm->isarr)) ) {
    if ( (err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, gtype,fdm->isarr)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from fd_setcompdom, err=%d\n",err);
	fflush(stream);
	return err;
    }
    /*--------------------------------------------------------------------------*/
    /*-declare computational domain---------------------------------------------*/
    err = rd_a_declare(&(model->ld_c), RDOM_MAX_NARR, ndim, dgs, dge);
    if ( err ) {
	fprintf(stream, "ERROR. fd_modelcrea from rd_a_declare err=%d\n", 
		err);
	fflush(stream);
	return E_BADINPUT;
    }
    /*--------------------------------------------------------------------------*/
    /*-compute send and receive domains------------------------------------*/
    for ( idim = 0; idim < RDOM_MAX_NARR; ++idim ) {
	IASN(dgsa[idim], IPNT_1);
	IASN(dgea[idim], IPNT_0);
    }
    for ( iv = 0; iv < RDOM_MAX_NARR; ++iv ) {
	//    fprintf(stderr,"iv=%d\n",iv);
	err = ex_compute(iv, &sten, &(model->ld_c), 
			 dgsa[iv], dgea[iv], dgsr, dger, frcvempty);
	if ( err ) {
	    fprintf(stream, "ERROR. fd_modelcrea from ex_compute err=%d for array [%s].\n", 
		    err, (fdm->ind2str)(iv));
	    fflush(stream);
	    return E_INTERNAL;
	}
	/* TODO: change process receives below */
	/* First, check that only P and V are received */
	rcvne = 0;
	for ( i = 0; i < nnei; ++i ) if ( !(frcvempty[i]) ) ++rcvne;
            
	/* check which array exchanged */
	if ( rcvne > 0 ) {
	    for ( idim = 0; idim < ndim; ++idim ) 
		if ( fdm->isdyn(iv) ) break;
	    if ( idim == ndim ) {
		fprintf(stream, "ERROR. fd_modelcrea: wrong array to be exchanged [%s].\n", 
			fdm->ind2str(iv));
		fflush(stream);
		return E_INTERNAL;
	    }
	}
	/* Second, store receive domains */
	for ( i = 0; i < nnei; ++i ) {
	    //      fprintf(stream,"iv=%d i=%d dgsr[i][0]=%d dger[i][0]=%d dgsr[i][1]=%d dger[i][1]=%d\n",
	    //	      iv, i, dgsr[i][0],dger[i][0],dgsr[i][1],dger[i][1]);
	    IASN(dgsrs[i][iv], dgsr[i]);
	    IASN(dgers[i][iv], dger[i]);
	}

	/* model-specific alterations of any arrays */
	//    fprintf(stderr,"fd_modelcrea->alter_dom,iv=%d\n",iv);
	err = fdm->alter_dom(iv,dgsa[iv],dgea[iv]);
	if (err) {
	    fprintf(stream, "ERROR. fd_modelcrea from alter_dom, err=%d\n", err);
	    fflush(stream);
	    return err;
	}      
	//    fprintf(stderr,"bottom iv=%d\n",iv);
    }
    /*--------------------------------------------------------------------------*/
    /*-allocate main domain, create computational domain------------------------*/
    err = rd_a_create(&(model->ld_a), RDOM_MAX_NARR, ndim, dgsa, dgea);
    if ( err ) {
	fprintf(stream, "ERROR. fd_modelcrea from rd_a_create allocated domain err=%d.\n", err);
	fflush(stream);
	return E_INTERNAL;
    }
    model->ld_c = model->ld_a;
    //  for (i = 0;i < fdm->getnarr();i ++) {
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    iv = fdm->getindices(i);
	//    if ( !(fdm->isdyn(iv)) ) continue;
	if ( !(fdm->isdyn(i)) ) continue;
	//    err = rd_greset(&(model->ld_c), iv, dgs[iv], dge[iv]);
	err = rd_greset(&(model->ld_c), i, dgs[i], dge[i]);
	if ( err ) {
	    fprintf(stream, 
		    "ERROR. fd_modelcrea from rd_greset computational array [%s] err=%d.\n", 
		    //	      fdm->ind2str(iv), err);
		    fdm->ind2str(i), err);
	    fflush(stream);
	    return E_INTERNAL;
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-create physical domain, here it includes the 2 bnd pnts------------------*/
    /*-not for sure, need advise from Dr. Symes!!!!!!---------------------------*/
    get_n(ns, model->g);
    model->ld_p = model->ld_c;
    //  for (i = 0;i < fdm->getnarr();i ++) {
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    iv = fdm->getindices(i);
	//    if ( (err=rd_gse(&(model->ld_p), iv, gs, ge)) ) {
	if ( (err=rd_gse(&(model->ld_p), i, gs, ge)) ) {
	    fprintf(stream, 
		    "ERROR. fd_modelcrea from rd_gse physical array [%s] err=%d.\n", 
		    //              fdm->ind2str(iv), err);
		    fdm->ind2str(i), err);
	    fflush(stream);
	    return E_INTERNAL;
	}
	for (idim = 0;idim < ndim;idim ++) {
	    //      if (gtype[iv][idim] == PRIMAL_GRID) {
	    if (gtype[i][idim] == PRIMAL_GRID) {
#if INCLUDE_BOUNDARY_PNTS
		/* include left and right boundary pnts */
		if ( gs[idim] < 0 )            gs[idim] = iwave_min(0, ge[idim] + 1);
		if ( ge[idim] > ns[idim] - 1 ) ge[idim] = ns[idim] - 1;
#else
		/* not include left and right boundary pnts */
		if ( gs[idim] < 1 )            gs[idim] = iwave_min(1, ge[idim] + 1);
		if ( ge[idim] > ns[idim] - 2 ) ge[idim] = ns[idim] - 2;
#endif     

		if ( ge[idim] < gs[idim] )     ge[idim] = gs[idim] - 1;
	    }
	    //      else if (gtype[iv][idim] == DUAL_GRID) {
	    else if (gtype[i][idim] == DUAL_GRID) {
		if ( gs[idim] < 0 )            gs[idim] = iwave_min(0, ge[idim] + 1); 
		if ( ge[idim] > ns[idim] - 2 ) ge[idim] = ns[idim] - 2;
		if ( ge[idim] < gs[idim] )     ge[idim] = gs[idim] - 1;
	    }
	    else {
		fprintf(stream, "ERROR. fd_modelcrea: undefined grid type: %d\n", 
			//                gtype[iv][idim]);
			gtype[i][idim]);
		fflush(stream);
		return E_INTERNAL;
	    }
	}
	//    if ( (err=rd_greset(&(model->ld_p), iv, gs, ge)) ) {
	if ( (err=rd_greset(&(model->ld_p), i, gs, ge)) ) {
	    fprintf(stream, 
		    "ERROR. fd_modelcrea from rd_greset physical array [%s] err=%d.\n", 
		    //              fdm->ind2str(iv), err);
		    fdm->ind2str(i), err);
	    fflush(stream);
	    return E_INTERNAL;
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-set virtual receive domains ---------------------------------------------*/
    for (inei = 0;inei < nnei;inei ++) {
	model->ld_r[inei] = model->ld_s[inei] = model->ld_a;
	//    for (i = 0;i < fdm->getnarr();i ++) {
	for (i = 0;i < RDOM_MAX_NARR;i ++) {
	    //      iv = fdm->getindices(i);
	    //      if (!(fdm->isdyn(iv)))  continue;
	    if (!(fdm->isdyn(i)))  continue;
	    //      err = rd_greset(model->ld_r + inei, iv, dgsrs[inei][iv], dgers[inei][iv]);
	    err = rd_greset(model->ld_r + inei, i, dgsrs[inei][i], dgers[inei][i]);
	    if ( err ) {
		fprintf(stream, 
			"ERROR. fd_modelcrea from rd_greset\n receive array [%s] if domain (%d) greset error #%d.\n iv=%d dgsrs[0]=%d dgers[0]=%d dgsrs[1]=%d dgers[1]=%d\n",
			//		fdm->ind2str(iv), inei, err,iv,dgsrs[inei][iv][0],dgers[inei][iv][0],dgsrs[inei][iv][1],dgers[inei][iv][1]);
			fdm->ind2str(i), inei, err,i,dgsrs[inei][i][0],dgers[inei][i][0],dgsrs[inei][i][1],dgers[inei][i][1]);
		fflush(stream);
		return E_INTERNAL;
	    }
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-set local grid ----------------------------------------------------------*/
    /* must have a variable defined on PRIMAL grid in every axis */
    //  for (i = 0;i < fdm->getnarr();i ++) {
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	//    iv = fdm->getindices(i);
	for (idim = 0;idim < ndim;idim ++) {
	    //      if (gtype[iv][idim] == DUAL_GRID)  break;
	    if (gtype[i][idim] == DUAL_GRID)  break;
	}
	if (idim == ndim){
	    //      err = rd_gse(&(model->ld_c), iv, gs, ge);
	    err = rd_gse(&(model->ld_c), i, gs, ge);

	    if ( err ) {
		fprintf(stream, "ERROR. fd_modelcrea from rd_gse allocated array [%s] err=%d.\n",
			//                fdm->ind2str(iv), err);
			fdm->ind2str(i), err);
		fflush(stream);
		return E_INTERNAL;
	    }
	    break;
	}
    }
    init_grid(&(model->gl),ndim);
    for ( idim = 0; idim < ndim; ++idim ) {
	model->gl.axes[idim].d = model->g.axes[idim].d;
	model->gl.axes[idim].o = model->g.axes[idim].o + model->gl.axes[idim].d * gs[idim];
	model->gl.axes[idim].n = ge[idim]-gs[idim]+1;
	model->gl.axes[idim].id = model->g.axes[idim].id;
    }
    fprintf(stream,"NOTE: local grid used to determine trace sampling:\n");
    fprint_grid(stream,model->gl);
    /*--------------------------------------------------------------------------*/

    /* deallocate stencil */
    sten_destroy(&sten);

    return 0;
}

int fd_modeldest(IMODEL * model) 
{
    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    return fdm->fd_model_dest(model);
}

int fd_isdyn(IMODEL * model, int i) {
    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    return fdm->isdyn(i);
}

int fd_modelinfo(FILE * stream, IMODEL * model) {
    return 0;
}

int fd_mread(FILE * stream, IMODEL * model, PARARRAY * par, int panelindex) {
    FD_MODEL *fdm = (FD_MODEL *)model->specs;
    int err;

    err = fdm->readmedia(par, stream, model, panelindex);
    if (err) {
	fprintf(stream, "Error: error in fd_mread after calling readmedia\n");
	return err;
    }

    err = fdm->readtimegrid(par, stream, model);
    if (err) {
	fprintf(stream, "Error: error in fd_mread after calling readtimegrid\n");
	return err;
    }
    if ((model->tsind).dt <= REAL_ZERO) {
	fprintf(stream, "Error: bad input: wrong time step dt=%g\n", 
		(model->tsind).dt);
	return E_BADINPUT;
    }
#ifdef VERBOSE
    fprintf(stderr, "dt = %g\n", (model->tsind).dt);
#endif
  
    //  fprintf(stderr,"mread->schemeinfo, ndim=%ld\n",(model->g).dim);
    //  err = fdm->readschemeinfo(par, stream, fdm, model);
    // changed interface 13.02.11 WWS
    err = fdm->readschemeinfo(par, stream, model);
    if (err) {
	fprintf(stream, "Error: error in fd_mread from readschemeinfo\n");
	return err;
    }
  
    return 0;
}

int im_destroy(IMODEL *model)
/*< Destroys model (STORAGE DEALLOCATION). >*/
{
  rd_a_destroy(&(model->ld_a));
  free(model->ld_s);
  fd_modeldest(model);
  free(model->specs);
  return im_construct(model);
}
/*----------------------------------------------------------------------------*/
