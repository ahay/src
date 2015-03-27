#ifndef __SEAM_GRID
#define __SEAM_GRID

/* tolerance parameter - in C++ would get
   from numeric limits. use: two grid locations
   x and y are identified if abs(x-y) < TOL*d.
*/
#define TOL 100.0*REAL_EPS

// define boundary index between external, internal extended axis indices
#define EXTINT 100

#include "utils.h"
#include "except.hh"

/** Axis: Basic Metadata element for regular grids. Defines a uniformly
    sampled axis, after the fashion of SEPlib77 or RSF. Struct part of
    RVLGrid Axis class. 
    
    Innovative feature: each axis carries an integer token (axis.id)
    meant to signify its position in a global ordering of axes (hence
    its physical meaning)..
*/

typedef struct s_axis {
  /** number of gridpoints on axis */
  size_t n; 
  /** step between gridpoints */
  ireal d;
  /** coordinate of first gridpoint */
  ireal o;
  /** axis index */
  int id;
} axis;

/** default constructor. Default values: n=1, d=1.0, o=0.0, id=1
 @param[out]  a (axis *) - axis to be initialized 
*/
void init_default_axis(axis * a);
  
/** main constructor - not terribly useful. id defaults to 1.
@param[out] a (axis *) - axis to be initialized 
@param[in] n (int) - number of gridpoints on axis
@param[in] d (ireal) - step
@param[in] o (ireal) - coord of first data point
*/
void init_axis(axis * a, size_t n, ireal d, ireal o);
  
/** deep copy function
@param[out] tgt - target axis
@param[in] src - source axis
*/
void copy_axis(axis * tgt, const axis * src);

/** print axis to stdout
@param[in] a (axis) - axis to be printed
 */
void print_axis(axis a);
  
/** print axis to stream 
@param[in] fp (FILE *) - output stream
f@param[in] a (axis) - axis to be printed
 */
void fprint_axis(FILE * fp, axis a);
void fprint_num_axis(FILE * fp, int i, axis a);
/** compare two axes for equality
@param[in] a1 (axis) - first axis 
@param[in] a2 (axis) - second axis
@return 0 if axes are same, else 1
*/
int compare_axis(axis a1, axis a2);
/** Regular grid struct 
@param dim - dimension of physical grid 
@param gdim - global dimension including nonphysical axes
@param axis - vector of axes defining grid
*/
typedef struct {
  int dim;
  int gdim;
  axis axes[RARR_MAX_NDIM];
} grid;
  
/** default constructor. 
@param[out] g (grid *) - grid to be initialized
*/
void init_default_grid(grid * g);

/** main constructor. initializes metadata (dim) but not data
    (axes). Since no dynamic allocation takes place, only real role is
    to sanity-check dim.
@param[out] g (grid *) - grid to be initialized
@param[in] dim (int) - dimension of grid, at most \ref RARR_MAX_NDIM
@return 0 for normal return, E_BADINPUT if input inconsistent
*/
int init_grid(grid * g, int dim, int gdim);

/** deep copy 
@param[out] tgt - target grid
@param[in] src - source grid
*/
void copy_grid(grid * tgt, const grid * src);

/** print grid to stdout 
@param[in] a (grid) - grid to be printed
*/
void print_grid(grid a);
  
/** print grid to stream 
@param[in] a (grid) - grid to be printed
@param[in] fp (FILE *) - stream to which to print
*/
void fprint_grid(FILE * fp, grid a);
  
/** compare two grids by comparing their axes 
@return 0 if grids same, else 1
*/
int compare_grid(const grid g1, const grid g2);

/** decide if two grids are compatible, that is, 
    define sublattices of a common parent lattice
    @return 0 if compatible, else 1
*/
int compatible_grid(const grid g1, const grid g2);

/** return effective (physical plus internal extended)
    grid dimension */
int get_dimension_grid(grid g);

/** get number of physical gridpoints (product of n's)
    @return product of physical axis lengths
 */
int get_datasize_grid(grid g);

/** get number of physical and internal extended gridpoints (product of n's)
    @return product of physical, internal extended axis lengths
 */
int get_extended_datasize_grid(grid g);

/** get total number of gridpoints (product of n's)
    @return product of axis lengths
 */
size_t get_global_datasize_grid(grid g);

/** returns (physical, non-extended) cell vol */  
ireal get_cellvol_grid(grid g);

/** returns physical + internal extended cell vol */
ireal get_extended_cellvol_grid(grid g);

/** returns global (all axes) cell vol */  
ireal get_global_cellvol_grid(grid g);

/** get total number of records = physical (+ internal ext'd) grids
    within global grid (product of n's) 
    @return product of external ext'd axis lengths
 */
int get_panelnum_grid(grid g);

/** get axis length array 
    @param[out] n (IPNT) - axis lengths
*/
void get_n(IPNT n, grid g);
/** get step array 
    @param[out] d (RPNT) - steps
*/
void get_d(_RPNT d, grid g);
/** get grid origin coordinate array 
    @param[out] o (RPNT) - grid origin
*/
void get_o(_RPNT o, grid g);

/** get array of indices of grid origin in global grid 
    @param[out] gs (IPNT) - global indices of grid origin
*/
void get_gs(IPNT gs, grid g);

/** get array of indices of grid origin in global grid 
    @param[out] gs (IPNT) - global indices of grid origin
*/
void get_ge(IPNT ge, grid g);

/** get array of axis ids as function of axis index - next function
    inverts this relation */
void get_id(IPNT id, grid g);

/** returns axis order array, i.e. axis index as a function of id,
    rather than id as a function of axis index (which is stored). 
    @param[out] a (IPNT) - axis order array
*/
void get_ord(IPNT a, grid g);

/** convex hull operation. Loops through axis vector input. If axis
    has new id, then added to grid. If axis id same as existing, then
    axis enlarged convex hull of existing and new axis. Axes must be 
    commensurable.
*/
bool grid_union(grid * g, axis const * ax);

/** initializes steps through the non-spatial axes of a
    grid. Assumption is that time is the first axis, so should be
    reversed for adjoint stepping.
 */
bool init_step(grid g, IPNT step, bool fwd);

/** increments index array for next step in simulation - time axis
    first, then other non-spatial axes in order. Time axis decremented
    for adjoint case. Returns true until no more steps are possible
    (including along non-time axes), then false. */
bool next_step(grid g, IPNT step);

#endif /* __SEAM_GRID */


