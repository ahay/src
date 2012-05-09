#ifndef __SEAM_GRID
#define __SEAM_GRID

/* tolerance parameter - in C++ would get
   from numeric limits. use: two grid locations
   x and y are identified if abs(x-y) < TOL*d.
*/
#define TOL 0.01
  
#include "utils.h"

/** Axis: Basic Metadata element for regular grids. Defines a uniformly
    sampled axis, after the fashion of SEPlib77 or RSF. Struct part of
    RVLGrid Axis class. 
    
    Innovative feature: each axis carries an integer token (axis.id)
    meant to signify its position in a global ordering of axes (hence
    its physical meaning)..
*/

typedef struct {
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
int init_default_axis(axis * a);
  
/** main constructor - not terribly useful. id defaults to 1.
@param[out] a (axis *) - axis to be initialized 
@param[in] n (int) - number of gridpoints on axis
@param[in] d (ireal) - step
@param[in] o (ireal) - coord of first data point
*/
int init_axis(axis * a, size_t n, ireal d, ireal o);
  
/** print axis to stdout
@param[in] a (axis) - axis to be printed
 */
int print_axis(axis a);
  
/** print axis to stream 
@param[in] fp (FILE *) - output stream
f@param[in] a (axis) - axis to be printed
 */
int fprint_axis(FILE * fp, axis a);
  
/** compare two axes for equality
@param[in] a1 (axis) - first axis 
@param[in] a2 (axis) - second axis
@return 0 if axes are same, else 1
*/
int compare_axis(axis a1, axis a2);
 
/** Regular grid struct */
typedef struct {
  size_t dim;
  axis axes[RARR_MAX_NDIM];
} grid;
  
/** default constructor. 
@param[out] g (grid *) - grid to be initialized
*/
int init_default_grid(grid * g);

/** main constructor. initializes metadata (dim) but not data
    (axes). Since no dynamic allocation takes place, only real role is
    to sanity-check dim.
@param[out] g (grid *) - grid to be initialized
@param[in] dim (int) - dimension of grid, at most \ref RARR_MAX_NDIM
*/
int init_grid(grid * g, size_t dim);

/** print grid to stdout 
@param[in] a (grid) - grid to be printed
*/
int print_grid(grid a);
  
/** print grid to stream 
@param[in] a (grid) - grid to be printed
@param[in] fp (FILE *) - stream to which to print
*/
int fprint_grid(FILE * fp, grid a);
  
/** compare two grids by comparing their axes 
@param[in] g1 (grid) - first grid
@param[in] g2 (grid) - second grid
@return 0 if grids same, else 1
*/
int compare_grid(grid g1, grid g2);
  
/** get number of gridpoints (product of n's)
@param[in] g (grid) - input grid
@return product of axis lengths
 */
int get_datasize_grid(grid g);
  
/** get axis length array 
@param[in] g (grid) - input grid
@param[out] n (IPNT) - axis lengths
*/
int get_n(IPNT n, grid g);
/** get step array 
@param[in] g (grid) - input grid
@param[out] d (RPNT) - steps
*/
int get_d(_RPNT d, grid g);
/** get grid origin coordinate array 
@param[in] g (grid) - input grid
@param[out] o (RPNT) - grid origin
*/
int get_o(_RPNT o, grid g);
/** get array of indices of grid origin in global grid 
@param[in] g (grid) - input grid
@param[out] gs (IPNT) - global indices of grid origin
*/
int get_gs(IPNT gs, grid g);

/** returns axis order array, i.e. axis index as a function of id,
    rather than id as a function of axis index (which is stored). 
@param[in] g (grid) - input 9grid
@param[out] a (IPNT) - axis order array
*/
int get_ord(IPNT a, grid g);

#endif /* __SEAM_GRID */


