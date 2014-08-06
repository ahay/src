#ifndef __IWAVE_MOVIE_OUTPUT__
#define __IWAVE_MOVIE_OUTPUT__

#define KEYLEN 16
#define DEFSTEP 100

#include "traceio.h"
#include "gridio.h"
#include "model.h"

/** \section movie Movie

IWAVE movie facility. Provides ability to construct 2D movies from either 2D or 3D 
simulation. Any fields in the RDOM defining the simulation may be sampled; codes
for these fields and their indices in the RDOM should be defined in a subclass 
for a particular model. Extracts 2D plane from 3D volumes. Allows choice of plane 
(parallel to a coordinate plane) and position of plane in volume.

*/

typedef struct {
  int frameindex;
  int framestep;
  int framestart;
  int framestop;
  int it;
  int dim3d;           /* slice axis (3D) */
  int slice3d;         /* slice index (3D) */
  int maxslice;        /* max size of slice - for buffer alloc */
  ireal * buf;         /* buffer for use in noncontig case */
  char* hname;         /* also used as flag to create movies */
  grid mg;             /* movie grid for output */
  int imovie[RDOM_MAX_NARR];     /* RDOM indices */
  char * smovie[RDOM_MAX_NARR];  /* filenames */
  int nmovie;                    /* number of movies */

  /** to be assigned by constructor of concrete child class - 
      @param[in] key (const char *) - key for field, legit values
      defined in child class. Example: for pressure-velocity acoustics,
      fields that might be sampled are p and the three velocity components.
      Accordingly, a child class could assign this pointer to a function
      which returns appropriate RDOM indices for keys p, v1, v2, v3.

      @return >=0 and <= RDOM_MAX_NARR indicates successful return, 
      gives index into rdom; < 0 indicates invalid key
  */
  int (*iselect)(const char * key);

} MOVIE;

/** extract a snapshot if time is one of specified snapshot times.
    @return success=0
 *  @param[in]   _irec (int) - current record index (now take 0 for standard modeling; Not support extended modeling) 
*/

int movie_run(MOVIE *t, IMODEL *m, FILE * stream, int _irec);

/** movie facility initialization controlled by PARARRAY (3rd arg). Parameters are:
    <ol>
    <li> key = movie[n], value =  string, signifying field to be sampled. Any 
    number up to the number of available fields may be listed. Example: movie1=p, 
    movie2=v1. </li>
    <li> key = moviestep, value = real, step in ms between movie frames. Default = 
    DEFSTEP * internal timestep. DEFSTEP defined in this file.  </li>
    <li> key = movieaxis3d, value = int, for 3D only: only 2d movies possible - 
    indicates axis perpindicular to movie plane. Default = 2, i.e. the plane of 
    the first two coordinates will be sampled</li>
    <li> key = movieslice3d, value = read, for 3D only: position along axis at 
    which to position movie plane (perp. to axis). Default = 0, i.e. boundary 
    face will be sampled.</li>
    </ol>
    All params are defaulted, except for movie[n] (i.e. movie1, movie2, movie3,...). 
    If movie[n] param not present, or only gives invalid field ids, then no movie 
    recorded. 

    @return = 0 indicates successful initialization
*/
int movie_init(MOVIE *t, IMODEL *m, PARARRAY *par, tracegeom *tg, FILE *stream);

/** prints out movie parameters */
void movie_fprint(MOVIE *t, FILE *fp);

/** destructor 

    @return success=0
*/
int movie_destroy(MOVIE *t);

/** default initializer: should be called first in any concrete
    subclass constructor.

    prototype for the constructor - MOVIE struct represents an 
    abstract base so doesn't have one:

    int mykindamovie_construct(MOVIE *t,  FILE *stream);

    a constructor should (a) call setnull, and (b) assign iselect to
    an actual function pointer.

    @return success=0
*/
int movie_setnull(MOVIE *t);
#endif
