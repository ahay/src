/* Functions for extraction of offset lists of subarrays.

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

#include <trip/base.h>

#include "offsets.h"

int get_array_offsets(off_t ** offs,
		      size_t * noffs,
		      int dim,
		      _IPNT gs,
		      _IPNT gn,
		      _IPNT ls,
		      _IPNT ln) 
/*< ** @fn get_array_offsets(off_t ** offs,
  size_t * noffs,
  int dim,
  IPNT gs,
  IPNT gn,
  IPNT ls,
  IPNT ln);
  @brief computes offsets into a subarray.
  
  same as @see get_array_offsets,
  
  but uses off_t, size_t instead of explicit types.
>*/
{

  /********************
   *   declarations   *
   ********************/
  
  int err=0;
  /* loop counters */
  size_t i,j,k;          
  /* dim-1 dim'l loc power sequence excluding axis 0 */
  size_t ploc[RARR_MAX_NDIM];
  /* dim dim'l glob power sequence */
  size_t pglob[RARR_MAX_NDIM];
  /* storage for expansion coefficients */
  int idx[RARR_MAX_NDIM];

  off_t a;
  off_t b;

  /********************
   * end declarations *
   ********************/

#ifdef DEBUG_OFFSETS
  for (i=0;i<RARR_MAX_NDIM;i++)
    fprintf(stderr,"i=%d gs[i]=%d gn[i]=%d ls[i]=%d ln[i]=%d\n",i,gs[i],gn[i],ls[i],ln[i]);
#endif

  /* bail if loc is empty */
  for (i=0; i<dim; i++)
    err = err || (ls[i]<gs[i]) || (ls[i]+ln[i]>gs[i]+gn[i]);
  if (err) {
#ifdef DEBUG_OFFSETS
    fprintf(stderr,"Error: get_array_offsets\n");
    fprintf(stderr,"local array is not subarray of global array\n");
#endif
    return E_OUTOFBOUNDS;
  }

  /* sanity check */
  if (dim > RARR_MAX_NDIM) {
#ifdef DEBUG_OFFSETS
    fprintf(stderr,"Error: get_array_offsets\n");
    fprintf(stderr,"dim = %d larger than permissible limit %d\n",
	    dim,RARR_MAX_NDIM);
#endif
    return E_OUTOFBOUNDS;
  }

  /* work out number of offsets */
  *noffs=1;
  for (i=1; i<dim ;i++) *noffs *= ln[i];

  /* it's convenient to accommodate the degenerate case noffs=0,
     caused by one of the ln's=0. Just return offs=NULL. */
  
  if (*noffs<1) {
    offs=NULL;
  }
  else {

    /* initialize workspace */
    *offs = (off_t *)malloc(*noffs*sizeof(off_t));

    /* compute power sequences. Note that ploc is dim-1 dim'l, and 
       consists of (1,ln1,ln1*ln2,...), whereas pglob is dim dim'l and
       consists of (1,n0,n0*n1,...).*/
    ploc[0]=1;
    pglob[0]=1;
    for (i=1;i<dim-1;i++) {
      ploc[i]=ploc[i-1]*ln[i];
      pglob[i]=pglob[i-1]*gn[i-1];
    }
    if (dim>1) pglob[dim-1]=pglob[dim-2]*gn[dim-2];

#ifdef DEBUG_OFFSETS
    for (i=0;i<dim-1;i++) fprintf(stderr,"offsets: ploc[%d]=%lu\n",i,ploc[i]);
    for (i=0;i<dim;i++) fprintf(stderr,"offsets: pglob[%d]=%lu\n",i,pglob[i]);
#endif

    /* now loop through the enumerated offsets. Expand each index
       in the local power sequence (store in idx), then use these
       coefficients to expand the global offset in the global 
       power sequence. Use the axis 0 loop offset explicitly to 
       initiate the expansion, as only axes 1 and up are represented
       in idx. */
    for (j=0;j<*noffs;j++) {
      k=j;
      for (i=dim-1;i>0;i--) {
	idx[i]=k/ploc[i-1];
	k=k-ploc[i-1]*idx[i];
#ifdef DEBUG_OFFSETS
	fprintf(stderr,"j=%jd i=%jd k=%jd idx=%d\n",j,i,k,idx[i]);
#endif
      }
      idx[0]=k; /* should = 0 , since ploc[0]=1 ex def */
      (*offs)[j]=0;
      for (i=0;i<dim;i++) {
	a=pglob[i];
	b=(idx[i]+ls[i]-gs[i]);
	(*offs)[j]+=pglob[i]*b;
	/*	(*offs)[j]+=(idx[i]+ls[i]-gs[i])*pglob[i];*/
#ifdef DEBUG_OFFSETS
	fprintf(stderr,"offset comp: i=%d ls=%d gs=%d pglob=%ld off_tmp=%jd\n",i,ls[i],gs[i],pglob[i],(*offs)[j]);
#endif
      }

#ifdef DEBUG_OFFSETS
      fprintf(stderr,"offset: offs[%d]=%jd\n",j,(*offs)[j]);
#endif
    }
  }
  
  return err;
}
