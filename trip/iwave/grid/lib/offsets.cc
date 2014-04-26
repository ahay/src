#include "offsets.h"

//#define DEBUG_OFFSETS

int get_array_offsets(off_t ** offs,
		      size_t * noffs,
		      int dim,
		      const _IPNT gs,
		      const _IPNT gn,
		      const _IPNT ls,
		      const _IPNT ln) {

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

  off_t b;

  /********************
   * end declarations *
   ********************/

#ifdef DEBUG_OFFSETS
  for (i=0;i<dim;i++)
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
    *offs = (off_t *)usermalloc_(*noffs*sizeof(off_t));

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
