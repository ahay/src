#ifndef __SEAM_OFFSETS
#define __SEAM_OFFSETS

#include "utils.h"

/*
#define DEBUG_OFFSETS
*/
/** @file
    @brief Functions for extraction of offset lists of subarrays.

   Compute offsets for a subarray with starts ls and lengths ln of a
   global array with starts gs and lengths gn. Dimension of both
   assumed to be dim.

   number of offsets = product of ln[1]*...*ln[dim-1]. 

   For each offset j,

   <UL>
   <LI> figure out the index array idx:

   idx[dim-1]=j/(ln[1]*...*ln[dim-2]), remainder
   k=j-idx[dim-1]*ln[1]*...*ln[dim-2]

   for j=dim-2,j>0: idx[j]=k/(ln[1]*...*ln[j-1]),
   k=k-idx[j]*ln[1]*...*ln[j-1]
   </LI>
   
   <LI> construct the offset in the global array:

   offs[j]=sum from i=0:dim-1
   (idx[i]+ls[i]-gs[i])*gn[0]*...*gn[i-1]
   </LI>
   </UL>

   Example: suppose dim=2, gs=[0,0], gn=[4,4], ls=[1,1],
   ln=[2,2]. The subarray consists of the four points
   [1,1],[2,1],[1,2], and [2,2]. The offsets are clearly 5 and 9,
   noffs=2 = ln[1]*...*ln[dim-1].

   To work out the 2nd offset (j=1), product of ln[dim-2]*...*ln[1]=1
   (empty products are by default = 1). 1/1=1, remainder is 0. So
   idx=[0,1].

   2nd offset = (idx[0]+ls[0]-gs[0]) 
   (idx[1]+ls[1]-gs[1])*gn[0]
   = (0+1-0)+(1+1-0)*4=9 
*/

/** @fn get_array_offsets(off_t ** offs,
                          size_t * noffs,
                          int dim,
			  IPNT gs,
			  IPNT gn,
			  IPNT ls,
			  IPNT ln);
    @brief computes offsets into a subarray.
    
    same as @see get_array_offsets,

    but uses off_t, size_t instead of explicit types.
*/

int get_array_offsets(off_t ** offs,
		      size_t * noffs,
		      int dim,
		      const IPNT gs,
		      const IPNT gn,
		      const IPNT ls,
		      const IPNT ln);

#endif
