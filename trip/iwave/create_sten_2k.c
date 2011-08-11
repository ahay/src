/* Create stencils of a 2-2k FD scheme */
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

#include <stdio.h>
/*^*/

#include <trip/base.h>
/*^*/

#include "stencil.h"
/*^*/

#include "create_sten_2k.h"
#include "fd.h"

int create_sten2_2k(FILE * stream, 
		    //		    int k, int ndim, int m_size,
		    int k, int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    int (*isdyn)(int), 
		    //		    int (*getindices)(int),
		    STENCIL * sten ) 
/*<*
 * Create stencils of a 2-2k FD scheme for a 1st order wave equations provided the 
 * dependent matrix and grid types of variables (static and dynamic) are given.
 *
 * @param [in] stream - (FILE *) file stream for dumping run-time messages and warnings
 * @param [in] k      - (int) spatial order 2k
 * @param [in] ndim   - (int) grid dimension
 * @param [in] m_size - (int) number of arrays in rdomain
 * @param [in] gtype  - (IPNT *) grid type array, length = m_size
 * @param [in] stendm - (int **) stencil dependency matrix, m_size x m_size
 * @param [in] isdyn  - (int (*)(int)) - fcn ptr: returns true if index pts to dyn arr
 * @param [out] sten  - (STENCIL *) STENCIL struct pointer
 >*/
{
  
  int err = 0;
  int idim, i, j, ir, ip, next, iv, nmask = 0;
  int len, disp;
  STENCIL_MASK mask;
  IPNT ind;
  int val;
  int sten_dep_type;

  sten_setnull(sten);
  if (k < 1) return E_BADINPUT;
  
  for (i = 0;i < RDOM_MAX_NARR;i ++) {
    //    ir = getindices(i);
    //    if ( !( isdyn(ir) ) ) continue;
    if (!(isdyn(i))) continue;
    for (j = 0;j < RDOM_MAX_NARR;j ++) {
      //      val = sten_dep_mat[ir][j];
      val = sten_dep_mat[i][j];
      while (val > 0) {
        nmask ++;
        val /= 10;
      }
    }
  }
  if (nmask == 0) {
    fprintf(stream, "Error: bad input in create_sten2_2k: nmask = %d\n",
            nmask);
    return E_BADINPUT;
  }
    
  if ( (err = sten_create(sten, nmask)) ) {
    return err;
  }

  next = 0;
  for (i = 0;i < RDOM_MAX_NARR;i ++) 
  {
    //    ir = getindices(i);
    //    if ( !isdyn(ir) ) continue;
    if ( !isdyn(i) ) continue;
    for (j = 0;j < RDOM_MAX_NARR;j ++) 
    {
      //      ip = getindices(j);
      //      val = sten_dep_mat[ir][ip];
      val = sten_dep_mat[i][j];
      while (val > 0) {
        /* extract each digit in sten_dep_type */
        sten_dep_type = val % 10;
        val = val / 10;

        IASN(ind, IPNT_0);
        idim = -1;
        
        if (sten_dep_type == DEP_F) {
          for (iv = 0;iv < RARR_MAX_NDIM;iv ++) {
	    //            if (gtype[ir][iv] != gtype[ip][iv]) {
            if (gtype[i][iv] != gtype[j][iv]) {
              fprintf(stream, "Error: bad input in fd_create_sten2_2k: ");
              fprintf(stream, "array %d is dependent of array %d, ", ir, ip);
              fprintf(stream, "but they are defined on different types of grids\n");
              sten_destroy(sten);
              return E_BADINPUT;
            }
          }
          len = 1;  disp = 0;  idim = 0;
        }
        else {
          if (sten_dep_type == DEP_DFDZ) {
            idim = 0;
            if (ndim < 1) idim = -1;
          }
          if (sten_dep_type == DEP_DFDX) {
            idim = 1;
            if (ndim < 2) idim = -1;
          }
          if (sten_dep_type == DEP_DFDY) {
            idim = 2;
            if (ndim < 3) idim = -1;
          }
          if (idim == -1) {
            fprintf(stream, "Error: bad input in fd_create_sten2_2k: in %d undefined sten type: %d \n", 
                    ndim, sten_dep_type);
            sten_destroy(sten);
            return E_BADINPUT;
          }
	  //          if (gtype[ir][idim] == gtype[ip][idim]) {
          if (gtype[i][idim] == gtype[j][idim]) {
            len = 2*k + 1;
            disp = 0;
          }
          else {
            len = 2*k;
	    //            if (gtype[ir][idim] == DUAL_GRID) disp = 1;
            if (gtype[i][idim] == DUAL_GRID) disp = 1;
            else disp = 0;
          }
          
        }
	//        if ( (err = mask_create(&mask, ip, ir, len)) ) {
        if ( (err = mask_create(&mask, j, i, len)) ) {
          sten_destroy(sten);
          return err;
        } 
        for (iv = 0;iv < len;iv ++) {
          ind[idim] = iv - (int)(len/2) + disp;
          if ( mask_set(&mask, iv, ind) ) {
            sten_destroy(sten);
	    return  E_INTERNAL;
          }
        }
        if ( sten_set(sten, next++, &mask) ) {
          sten_destroy(sten);
          return  E_INTERNAL;
        }
      }
    }
  }

  return 0;
}
