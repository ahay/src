#include "create_sten_2k.h"

int create_sten2_2k(FILE * stream,
		    IWaveInfo const & ic,
		    int k, int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten ) {
  
  int err = 0;
  int idim, i, j, next, iv, nmask = 0;
  int len, disp;
  STENCIL_MASK mask;
  IPNT ind;
  int val;
  int sten_dep_type;

  sten_setnull(sten);
  if (k < 1) return E_BADINPUT;
  
  for (i = 0;i < RDOM_MAX_NARR;i ++) {
    /*    ir = getindices(i); */
    /*  if ( !( isdyn(ir) ) ) continue; */
    //	if (!(isdyn(fdm,i))) continue;
    if (!(fd_isdyn(i,ic))) continue;
    for (j = 0;j < RDOM_MAX_NARR;j ++) {
      /*      val = sten_dep_mat[ir][j]; */
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
      //	if ( !isdyn(fdm,i) ) continue;
      if ( !fd_isdyn(i,ic) ) continue;
      for (j = 0;j < RDOM_MAX_NARR;j ++) 
	{
	  val = sten_dep_mat[i][j];
	  while (val > 0) {
	    /* extract each digit in sten_dep_type */
	    sten_dep_type = val % 10;
	    val = val / 10;

	    IASN(ind, IPNT_0);
	    idim = -1;
        
	    if (sten_dep_type == DEP_F) {
	      for (iv = 0;iv < RARR_MAX_NDIM;iv ++) {
		/*            if (gtype[ir][iv] != gtype[ip][iv]) { */
		if (gtype[i][iv] != gtype[j][iv]) {
		  fprintf(stream, "Error: bad input in fd_create_sten2_2k: ");
		  /*              fprintf(stream, "array %d is dependent of array %d, ", ir, ip); */
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
	      /*          if (gtype[ir][idim] == gtype[ip][idim]) { */
	      if (gtype[i][idim] == gtype[j][idim]) {
		len = 2*k + 1;
		disp = 0;
	      }
	      else {
		len = 2*k;
		/*            if (gtype[ir][idim] == DUAL_GRID) disp = 1; */
		if (gtype[i][idim] == DUAL_GRID) disp = 1;
		else disp = 0;
	      }
          
	    }
	    /*       if ( (err = mask_create(&mask, ip, ir, len)) ) { */
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
