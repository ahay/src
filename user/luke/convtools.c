/* tools for convolutions */
/*^*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int conv_in_bounds( int *Ind, int *N, int ndim )
	/*< determines if an index is in bounds, return 0 if in bounds, 1 if out of bounds >*/
{
	/* looping index */
	int i;
	/* loop through to see if we go out of bounds */
	for ( i = 0 ; i < ndim ; i++ ){
		if ( Ind[ i] < 0 ) return 1;
		if ( Ind[ i] > N[ i] - 2 ) return 1;
	}
	return 0;
}

int conv_int_exponent(int j, int n)
	/*< returns j^n for integers j, n >= 0 >*/
{
	/* check to see if zeroth power */
	if ( n == 0 ) return 1;
	/* intialize */
	int out = 1;
	/* counting index */
	int i;
	/* loop through */
	for ( i = 0 ; i < n ; i++ ){
		/* multiply by j again */
		out *= j ;
	}
	/* return value */
	return out;
}

float* conv_scale_float_array( float *arrayin, float scale, int n)
	/*< scales ar array of size n >*/
{
	/* declare array */
	float* arrayout = sf_floatalloc(n);
	/* looping index */
	int i;
	for ( i = 0 ; i < n ; i++){
		arrayout[ i] = scale * arrayin[ i];
	}
	return arrayout;
}

float* conv_float_array_subtract( float *A, float *B, int ndim)
	/*< outputs C = A - B for float arrays A, B >*/
{
	/* declare C */
	float* C = sf_floatalloc(ndim);
	/* looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++ ){
		C[ i] = A[ i] - B[ i];
	}
	return C;
}

float* conv_float_array_add( float *A, float *B, int ndim)
	/*< outputs C = A + B for float arrays A, B >*/
{
	/* declare C */
	float* C = sf_floatalloc(ndim);
	/* looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++ ){
		C[ i] = A[ i] + B[ i];
	}
	return C;
}

int* conv_int_array_add( int *A, int *B, int ndim)
	/*< outputs C = A + B for int arrays A, B >*/
{
	/* declare C */
	int* C = sf_intalloc(ndim);
	/* looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++ ){
		C[ i] = A[ i] + B[ i];
	}
	return C;
}

long conv_arraysize(int *N, int ndim)
	/*< returns the size of the array by unwrapping N >*/
{
	/* initialize */
	long nelements = 1;
	int dim;
	/* loop thru dimensions */
	for (dim = 0 ; dim < ndim ; dim++){
		/* increase total element size */
		nelements *= (long)N[dim];
	}
	return nelements;
}


int* conv_unwrap(long indx, int *N, int ndim)
	/*< unwraps index into coordinates >*/
{
    int* Ind = sf_intalloc(ndim);
	int j, n;
	long r = indx;
	for (j = ndim-1 ; j >= 0 ; j--){
        n = conv_arraysize(N,j);
		Ind [ j] = (int)((double)r / (double)n);
		r -= (long) Ind[j]*n;
	}
	return Ind;
}

long conv_wrap(int *Ind, int *N, int ndim)
/*< wraps input coordinates back into index >*/
{
	long indx = 0;
	int j, n; 
	for ( j = 0; j < ndim ; j++){
        n = conv_arraysize(N,j);
		indx += n * Ind[j];
	}
	return indx;
}

float* conv_index_to_coordinates(int *Ind, float *D, float *O, int ndim)
	/*< transforms an index location to a coordinate >*/
{
	/* allocate coordinate array */
	float* Coord = sf_floatalloc(ndim);
	/* declare looping index */
	int i ;
	for (i = 0 ; i < ndim ; i++){
		Coord[ i] = ((float) Ind[ i]) * D[ i] + O[ i];
	}
	return Coord;
}



int* conv_coordinates_to_index( float *Coord, float *D, float *O, int ndim)
	/*< converts coordinate array to nearest index rounding down >*/
{
	/* declare index array */
	int* Ind = sf_intalloc(ndim);
	/* looping index */
	int i ;
	/* loop through dimensions */
	for ( i = 0 ; i < ndim ; i++){
		Ind[ i] = (int)((Coord[ i] - O[ i])/D[i]);
	}
	return Ind;
}

float* conv_index_coords_remainder( int *Ind, float *Coord, float *D, float *O, int ndim)
	/*< finds the interpolation proportion remainder between the index and the coordinate >*/
{
	/* declare remainder array */
	float* Rem = sf_floatalloc(ndim);
	/* declare converted coord array */
    float* CoordC = sf_floatalloc(ndim);
	/* get converted coordinates */
	free (CoordC);
	CoordC = conv_index_to_coordinates( Ind, D, O, ndim);
	/* declare looping index */
	int i ; 
    /* find coordinate for index, subtract from coordinate */
	float* Diff = conv_float_array_subtract( Coord, CoordC, ndim);
	/* divide the difference by the increment */
	for ( i = 0 ; i < ndim ; i++){
		Rem [ i] = Diff[ i]/D[ i];
	}
	free (Diff);
	free (CoordC);
	return Rem ;
}



int* conv_Nnode(int ndim)
	/*< returns a N array for the nodes where each n is 2 >*/
{
	/* declare Nnode array */
	int* Nnode = sf_intalloc(ndim);
	/* declare looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++){
		Nnode[ i] = 2;
	}
	return Nnode;
}

int* conv_doughnut_wrap( int *Ind, int *N, int ndim)
	/*< checks to see if index is out of bounds.  If it is, performs doughnut wrapping >*/
{
	/* wrapped index */
	int* IndW = sf_intalloc(ndim);
	/* looping index */
	int i ;
	/* how many wraps off are we?*/
	int wraps;
	/* loop through dimensions */
	for ( i = 0 ; i < ndim ; i++ ){
		if ( Ind[ i] < 0 ){ 
			/* too low */ 
			/* how far too low */
			wraps = -1*Ind[ i] / N[ i] + 1;
			/* correct */
			IndW[ i] = Ind[ i] +  wraps * N[ i] ;

		}else{
			if ( Ind[ i] > N[ i] - 1){
				/* too high */
				/* how far too low */
				wraps = Ind[ i] / N[ i] ;
				/* correct */
				IndW[ i] = Ind[ i] -  wraps * N[ i] ;
			}else{
				/* Goldilox */
				IndW[ i] = Ind[ i];
			}
		}
	}
	
	return IndW ;
}

float conv_interpolation_weights( int *NInd, float *Rem, int ndim)
	/*< returns interpolation weights for Node Index and Remainder >*/
{
	/* intialize the weight */
	float weight = 1;
	/* dimensional looping index */
	int i ;
	for ( i = 0 ; i < ndim ; i++ ){
		/* check to see if we are "near" or "far" */
		if ( NInd[ i] == 0 ){ 
			/* near */
			weight *= ( 1. - Rem[ i] ); 
		} else {
			/* far */
			weight *= Rem[ i];
		}
	}
	return weight;
}

float conv_array_doughnut_interpolator( int *Ind1, float *Rem, float *array, int *N, int ndim )
	/*< interpolates between elements of array >*/
{
	/* value to be interpolated */
	float interp = 0;
	/* create N array for nodes */
	int* Nnode = conv_Nnode( ndim);
	/* Index array for node position */
	int* NInd = sf_intalloc( ndim);
	/* offset index array */
	int* Ind2 = sf_intalloc( ndim);
	/* and the in bounds index after doughnut wrapping */
	int* Ind3 = sf_intalloc( ndim);
	/* how many nodes are we dealing with ? */
	long nnodes = conv_arraysize( Nnode, ndim);
	/* looping index for nodes, using long for consistency with unwrapping programs */
	long nindx ;
	/* index for array position */
	long indx ;
	/* interpolation weights */
	float weight;
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		free (NInd);
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* get interpolation weight */
		weight = conv_interpolation_weights( NInd, Rem, ndim);
		/* continue if zero */
		if ( weight == 0. ) continue;
		/* add that to the current index */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* check to see if inbounds, if not, doughnut wrap */
		free (Ind3);
		Ind3 = conv_doughnut_wrap( Ind2, N, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind3, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		interp += array[ indx] * weight;
	}
	/* free unneeded arrays */
	free (Nnode);
	free (NInd);
	free (Ind2);
	free (Ind3);
	
	return interp;
}

float conv_array_interpolator( int *Ind1, float *Rem, float *array, int *N, int ndim )
	/*< interpolates between elements of array >*/
{
	/* value to be interpolated */
	float interp = 0;
	/* create N array for nodes */
	int* Nnode = conv_Nnode( ndim);
	/* Index array for node position */
	int* NInd = sf_intalloc( ndim);
	/* offset index array */
	int* Ind2 = sf_intalloc( ndim);
	/* how many nodes are we dealing with ? */
	long nnodes = conv_arraysize( Nnode, ndim);
	/* looping index for nodes, using long for consistency with unwrapping programs */
	long nindx ;
	/* index for array position */
	long indx ;
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		free (NInd);
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* add that to the current index */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind2, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		interp += array[ indx] * conv_interpolation_weights( NInd, Rem, ndim);
	}
	/* free unneeded arrays */
	free (Nnode);
	free (NInd);
	free (Ind2);
	return interp;
}

float* conv_array_adj_interpolator( int *Ind1, float *Rem, float interp, float *arrayout, int *N, int ndim )
	/*< adjoint interpolates between elements of array >*/
{
	/* create N array for nodes */
	int* Nnode = conv_Nnode( ndim);
	/* Index array for node position */
	int* NInd = sf_intalloc( ndim);
	/* offset index array */
	int* Ind2 = sf_intalloc( ndim);
	/* doughnut wrapped index */
	int* Ind3 = sf_intalloc( ndim);
	/* how many nodes are we dealing with ? */
	long nnodes = conv_arraysize( Nnode, ndim);
	/* looping index for nodes, using long for consistency with unwrapping programs */
	long nindx ;
	/* index for array position */
	long indx ;
	/* interpolation weights */
	float weight;
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		free ( NInd);
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* get interpolation weight, if zero continue */
		weight = conv_interpolation_weights( NInd, Rem, ndim);
		if ( weight == 0. ) continue;
		/* add that to the current index */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* doughnut wrap */
		free (Ind3);
		Ind3 = conv_doughnut_wrap( Ind2, N, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind3, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		arrayout[ indx] += interp * weight;
	}
	/* free unneeded arrays */
	free (Nnode);
	free (NInd);
	free (Ind2);
	free (Ind3);
	
	return arrayout;
}

float* conv_get_translations(long indx, float *trans, int *N, int ndim)
	/*< get local translation array >*/
{
	/* local translation array */
	float* X = sf_floatalloc(ndim);
	/* looping index */
	int i;
	/* how many elements in domain */
	long nelements = conv_arraysize(N,ndim);
	/* loop through */
	for ( i = 0 ; i < ndim ; i++ ){
		/* get coordinate */
		X[ i] = trans[ ((long)i) * nelements + indx];
	}
	return X;
}

float* conv_translate(float *arrayin, float *X, int *N, float *D, float *O, int ndim)
	/*< returns a translated version of the array by vector X.  this doesn't wrap >*/	
{
	/* looping index */
	long indx ; 
	/* position index */
	int* Ind1 = sf_intalloc(ndim);
	/* translated position index */
	int* Ind2 = sf_intalloc(ndim);
	/* determine translation array */
	int* TInd = conv_coordinates_to_index( conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* and remainder of translation for interpolation */
	float* TRem = conv_index_coords_remainder( TInd, conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize(N,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* loop through output array */
	for (indx = 0 ; indx < nelements ; indx++){
		/* determine where we are */
		free (Ind1);
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* check to see if we are in bounds */
		if ( conv_in_bounds( Ind2, N, ndim ) > 0 ) continue ; 
		/* interpolate, dont need += because looping through index */
		arrayout[ indx] = conv_array_interpolator( Ind2, TRem, arrayin, N, ndim );
	}
	/* free unneeded arrays */
	free (Ind1);
	free (Ind2);
	free (TInd);
	free (TRem);
	
	/* return translated array */
	return arrayout;
}

float* conv_translate_wrap(float *arrayin, float *X, int *N, float *D, float *O, int ndim, bool adj)
	/*< returns a translated version of the array by vector X, with doughnut wrapping so an adjoint >*/	
{
	/* looping index */
	long indx ; 
	/* position index */
	int* Ind1 = sf_intalloc(ndim);
	/* translated position index */
	int* Ind2 = sf_intalloc(ndim);
	/* determine translation array */
	int* TInd = conv_coordinates_to_index( conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* and remainder of translation for interpolation */
	float* TRem = conv_index_coords_remainder( TInd, conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize(N,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* loop through output array */
	for (indx = 0 ; indx < nelements ; indx++){
		/* determine where we are */
		free (Ind1);
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* write to output array */
		if (!adj){
			/* as interpolation */
			arrayout[ indx] = conv_array_doughnut_interpolator( Ind2, TRem, arrayin, N, ndim );
		}else{
			/* as adjoint interpoloation */
			arrayout = conv_array_adj_interpolator( Ind2, TRem, arrayin[ indx], arrayout, N, ndim );
		} 
	}
	/* free unneeded arrays */
	free (Ind1);
	free (Ind2);
	free (TInd);
	free (TRem);
	
	/* return translated array */
	return arrayout;
}

float* conv_var_translate_wrap(float *arrayin, float *trans, int *N, float *D, float *O, int ndim, bool adj)
	/*< returns a translated version of the array by variable translation, with doughnut wrapping so an adjoint >*/	
{
	/* looping index */
	long indx ; 
	/* position index */
	int* Ind1 = sf_intalloc(ndim);
	/* translated position index */
	int* Ind2 = sf_intalloc(ndim);
	/* translation array precursor */
	float* TIndpre = sf_floatalloc(ndim);
	/*  translation array */
	int* TInd = sf_intalloc(ndim);
	/* and translation remainder */
	float* TRem = sf_floatalloc(ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize(N,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* local translation array */
	float* X = sf_floatalloc(ndim);
	/* loop through input array, this is different from the constant case */
	for (indx = 0 ; indx < nelements ; indx++){
		/* read translation */
		free (X);
		X = conv_get_translations(indx, trans, N, ndim);	
		/* get translation precursor */
		free (TIndpre);
		TIndpre = conv_float_array_subtract( O, X, ndim);
		/* convert to translation index */
		free (TInd);
		TInd = conv_coordinates_to_index( TIndpre, D, O, ndim);
		/* get remainder */
		free (TRem);
		TRem = conv_index_coords_remainder( TInd, TIndpre, D, O, ndim);
		/* determine where we are */
		free (Ind1);
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		free (Ind2);
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* write to output array */
		if ( !adj ){
			/* as interpolation */
			arrayout[ indx] = conv_array_doughnut_interpolator( Ind2, TRem, arrayin, N, ndim );
		}else{
			/* as adjoint interpoloation */
			arrayout = conv_array_adj_interpolator( Ind2, TRem, arrayin[ indx], arrayout, N, ndim );
		} 
	}
	/* free unneeded arrays */
	free (Ind1);
	free (Ind2);
	free (TInd);
	free (TRem);
	free (X);
	/* return translated array */
	return arrayout;
}

int* conv_ker_shift(int *Kin, int *Nk, int ndim )
	/*< shifts the kernel index so it is centered about zero >*/
{
	/* allocate output kernel */
	int* Kout = sf_intalloc(ndim);
	/* dimension index */
	int i ;
	for ( i = 0 ; i < ndim ; i++){
		Kout [ i] = Kin [ i] - ( (Nk [ i] - 1)/2 );
	}
	return Kout;
}


float* conv_convolve_ker(float *arrayin, int *N, float *kernel, int *Nk, int ndim, bool adj)
	/*< convolves the input array with a kernel with odd number of elements in each dimension.  Assumes sampling in the two is the same, as well as their dimensionality >*/
{
	/* determine number of elements in array */
	long nelements = conv_arraysize( N,ndim);
	/* determine number of elements in kernel */
	long kelements = conv_arraysize(Nk,ndim);
	/* allocate output array */
	float* arrayout = sf_floatalloc(nelements);
	/* index for looping through outptut array */
	long indxA;
	/* index for looping through kernel */
	long indxK;
	/* coordinates corresponding to index on output array*/
	int* AInd = sf_intalloc(ndim);
	/* index for kernel before shifting */
	int* KIndpre = sf_intalloc(ndim);
	/* index corresponding to position on kernel */
	int* KInd = sf_intalloc(ndim);
	/* array index shifted by kernel */
	int* AKInd = sf_intalloc(ndim);
	/* we are doing integer shifts, so TRem is always zero*/
	float* TRem = sf_floatalloc(ndim);
	/* set to zero */
	TRem = conv_scale_float_array( TRem, 0. , ndim);
	/* loop through output array */
	for( indxA = 0 ; indxA < nelements ; indxA++){
		/* get position index*/
		free (AInd);
		AInd = conv_unwrap( indxA, N, ndim);
		/* loop through kernel */
		for ( indxK = 0 ; indxK < kelements ; indxK++ ){
			/* check to see if kernel nonzero */
			if ( kernel[ indxK] == 0 ) continue ;
			/* determine position in Kernel */
			free (KIndpre);
			KIndpre = conv_unwrap( indxK, Nk, ndim);
			/* shift to center */
			free (KInd);
			KInd = conv_ker_shift( KIndpre, Nk, ndim );
			/* shift Array Index by Kernel Position */
			free (AKInd);
			AKInd = conv_int_array_add( AInd, KInd, ndim);
			if ( !adj ){
				/* as interpolation */
				arrayout[ indxA] += kernel[ indxK] * conv_array_doughnut_interpolator( AKInd, TRem, arrayin, N, ndim );
			}else{
				/* as adjoint interpoloation */
				arrayout = conv_array_adj_interpolator( AKInd, TRem, arrayin[ indxA]*kernel[ indxK], arrayout, N, ndim );
			} 			
		}	
	}
	/* free unneeded arrays */
	free (AInd);
	free (KInd);
	free (AKInd);
	free (TRem);
	/* return finished product */
	return arrayout;
}

float* conv_convolve_ker_translate( float *arrayin, float *X, int *N, float *D, float *O, float *kernel, int *Nk, int ndim, bool adj)
	/*< translates by a fixed amount and convolves with a kernel >*/
{
	/* looping index in array */
	long indxA ; 
	/* looping index in kernel */
	long indxK ;
	/* position index */
	int* AInd1 = sf_intalloc(ndim);
	/* translated position index */
	int* AInd2 = sf_intalloc(ndim);
	/* determine translation array */
	int* TInd = conv_coordinates_to_index( conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* and remainder of translation for interpolation */
	float* TRem = conv_index_coords_remainder( TInd, conv_float_array_subtract( O, X, ndim), D, O, ndim);
	/* index in kernel prior to shifting */
	int* KIndpre;
	/* index corresponding to position on kernel */
	int* KInd = sf_intalloc(ndim);
	/* array index shifted by kernel */
	int* AKInd = sf_intalloc(ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize(N,ndim);
	/* determine number of elements in kernel */
	long kelements = conv_arraysize(Nk,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* loop through output array */
	for( indxA = 0 ; indxA < nelements ; indxA++){
		/* get position index*/
		free (AInd1);
		AInd1 = conv_unwrap( indxA, N, ndim);
		/* Translate by X */
		free (AInd2);
		AInd2 = conv_int_array_add( AInd1, TInd, ndim);
		/* loop through kernel */
		for ( indxK = 0 ; indxK < kelements ; indxK++ ){
			/* check to see if kernel nonzero */
			if ( kernel[ indxK] == 0 ) continue ;
			/* determine position in Kernel */
			free (KIndpre);
			KIndpre = conv_unwrap( indxK, Nk, ndim);
			/* shift to center */
			free (KInd);
			KInd = conv_ker_shift( KIndpre, Nk, ndim );
			/* shift Array Index by Kernel Position */
			free (AKInd);
			AKInd = conv_int_array_add( AInd2, KInd, ndim);
			if ( !adj ){
				/* as interpolation */
				arrayout[ indxA] += kernel[ indxK] * conv_array_doughnut_interpolator( AKInd, TRem, arrayin, N, ndim );
			}else{
				/* as adjoint interpoloation */
				arrayout = conv_array_adj_interpolator( AKInd, TRem, arrayin[ indxA]*kernel[ indxK], arrayout, N, ndim );
			} 			
		}	
	}	
	/* free unneeded arrays */
	free (AInd1);
	free (AInd2);
	free ( TInd);
	free ( TRem);
	free ( KInd);
	free (AKInd);
	
	return arrayout;
}


float* conv_convolve_ker_var_translate( float *arrayin, float *trans, int *N, float *D, float *O, 
    float *kernel, int *Nk, int ndim, bool adj)
	/*< translates by a variable amount and convolves with a kernel >*/
{
	/* looping index through output array*/
	long indxA ;
	/* looping index through kernel */
	long indxK ;  
	/* position index */
	int* AInd1 = sf_intalloc(ndim);
	/* translated position index */
	int* AInd2 = sf_intalloc(ndim);
	/* translation array precursor */
	float* TIndpre = sf_floatalloc(ndim);
	/*  translation array */
	int* TInd = sf_intalloc(ndim);
	/* and translation remainder */
	float* TRem = sf_floatalloc(ndim);
	/* kernel index prior to shifting */
	int* KIndpre = sf_intalloc(ndim);
	/* kernel index */
	int* KInd = sf_intalloc(ndim);
	/* array index shifted by kernel */
	int* AKInd = sf_intalloc(ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize( N,ndim);
	/* and in the kernel */
	long kelements = conv_arraysize(Nk,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* local translation array */
	float* X = sf_floatalloc(ndim);
	/* loop through output array */
	for( indxA = 0 ; indxA < nelements ; indxA++){
		/* read translation */
		free (X);
		X = conv_get_translations(indxA, trans, N, ndim);	
		/* get precursor for translation index */
		free (TIndpre);
		TIndpre = conv_float_array_subtract( O, X, ndim);
		/* convert to translation index */
		free (TInd);
		TInd = conv_coordinates_to_index( TIndpre, D, O, ndim);
		/* get remainder */
		free (TRem);
		TRem = conv_index_coords_remainder( TInd, TIndpre, D, O, ndim);		
		/* get position index*/
		free (AInd1);
		AInd1 = conv_unwrap( indxA, N, ndim);
		/* Translate by X */
		free (AInd2);
		AInd2 = conv_int_array_add( AInd1, TInd, ndim);
		/* loop through kernel */
		for ( indxK = 0 ; indxK < kelements ; indxK++ ){
			/* check to see if kernel nonzero */
			if ( kernel[ indxK] == 0 ) continue ;
			/* determine position in Kernel */
			free (KIndpre);
			KIndpre = conv_unwrap( indxK, Nk, ndim);
			/* shift to center */
			free (KInd);
			KInd = conv_ker_shift( KIndpre, Nk, ndim );
			/* shift Array Index by Kernel Position */
			free (AKInd);
			AKInd = conv_int_array_add( AInd2, KInd, ndim);
			if ( !adj ){
				/* as interpolation */
				arrayout[ indxA] += kernel[ indxK] * conv_array_doughnut_interpolator( AKInd, TRem, arrayin, N, ndim );
			}else{
				/* as adjoint interpoloation */
				arrayout = conv_array_adj_interpolator( AKInd, TRem, arrayin[ indxA]*kernel[ indxK], arrayout, N, ndim );
			} 			
		}	
	}		
	/* free unneeded arrays */
	free (AInd1);
	free (AInd2);
	free ( TInd);
	free ( TRem);
	free ( KInd);
	free (   X );
	free (AKInd);
	
	return arrayout;
}


float* conv_convolve_ker_var_translate_omp( float *arrayin, float *trans, int *No, float *D, float *O, 
    float *kernel, int *Nk, int ndim, bool adj)
	/*< translates by a variable amount and convolves with a kernel >*/
{
	/* looping index through output array*/
	long indxA ;
	/* looping index through kernel */
	long indxK ;  
	/* position index */
	int* AInd1 = sf_intalloc(ndim);
	/* translated position index */
	int* AInd2 = sf_intalloc(ndim);
	/* translation array precursor */
	float* TIndpre = sf_floatalloc(ndim);
	/*  translation array */
	int* TInd = sf_intalloc(ndim);
	/* and translation remainder */
	float* TRem = sf_floatalloc(ndim);
	/* kernel index prior to shifting */
	int* KIndpre = sf_intalloc(ndim);
	/* kernel index */
	int* KInd = sf_intalloc(ndim);
	/* array index shifted by kernel */
	int* AKInd = sf_intalloc(ndim);
	/* determine number of elements in array */
	long nelements = conv_arraysize( N,ndim);
	/* and in the kernel */
	long kelements = conv_arraysize(Nk,ndim);
	/* initialize output array */
	float* arrayout = sf_floatalloc(nelements);
	/* local translation array */
	float* X = sf_floatalloc(ndim);
	/* loop through output array */
#ifdef _OPENMP
#pragma omp parallel for private(indxA, X, TIndpre, TInd, TRem, AInd1, AInd2, indxK, KIndpre, KInd, AKInd)
#endif
	for( indxA = 0 ; indxA < nelements ; indxA++){
		/* read translation */
		free (X);
		X = conv_get_translations(indxA, trans, N, ndim);	
		/* get precursor for translation index */
		free (TIndpre);
		TIndpre = conv_float_array_subtract( O, X, ndim);
		/* convert to translation index */
		free (TInd);
		TInd = conv_coordinates_to_index( TIndpre, D, O, ndim);
		/* get remainder */
		free (TRem);
		TRem = conv_index_coords_remainder( TInd, TIndpre, D, O, ndim);		
		/* get position index*/
		free (AInd1);
		AInd1 = conv_unwrap( indxA, N, ndim);
		/* Translate by X */
		free (AInd2);
		AInd2 = conv_int_array_add( AInd1, TInd, ndim);
		/* loop through kernel */
		for ( indxK = 0 ; indxK < kelements ; indxK++ ){
			/* check to see if kernel nonzero */
			if ( kernel[ indxK] == 0 ) continue ;
			/* determine position in Kernel */
			free (KIndpre);
			KIndpre = conv_unwrap( indxK, Nk, ndim);
			/* shift to center */
			free (KInd);
			KInd = conv_ker_shift( KIndpre, Nk, ndim );
			/* shift Array Index by Kernel Position */
			free (AKInd);
			AKInd = conv_int_array_add( AInd2, KInd, ndim);


			if ( !adj ){
				/* as interpolation */
				arrayout[ indxA] += kernel[ indxK] * conv_array_doughnut_interpolator( AKInd, TRem, arrayin, N, ndim );
			}else{
				/* as adjoint interpoloation */
#ifdef _OPENMP
#pragma omp critical
#endif
				arrayout = conv_array_adj_interpolator( AKInd, TRem, arrayin[ indxA]*kernel[ indxK], arrayout, N, ndim );
			} 			
		}	
	}		
	/* free unneeded arrays */
	free (AInd1);
	free (AInd2);
	free ( TInd);
	free ( TRem);
	free ( KInd);
	free (   X );
	free (AKInd);
	
	return arrayout;
}