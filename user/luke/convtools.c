/* tools for convolutions */
/*^*/
#include <rsf.h>

int conv_in_bounds( int* Ind, int* N, int ndim )
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

float* conv_scale_float_array( float* arrayin, float scale, int n)
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

float* conv_float_array_subtract( float* A, float* B, int ndim)
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

float* conv_float_array_add( float* A, float* B, int ndim)
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

int* conv_int_array_add( int* A, int* B, int ndim)
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

long conv_arraysize(int* N, int ndim)
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


int* conv_unwrap(long indx, int* N, int ndim)
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

long conv_wrap(int* Ind, int* N, int ndim)
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

float* conv_index_to_coordinates(int* Ind, float* D, float* O, int ndim)
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



int* conv_coordinates_to_index( float* Coord, float* D, float* O, int ndim)
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

float* conv_index_coords_remainder( int* Ind, float* Coord, float* D, float* O, int ndim)
	/*< finds the interpolation proportion remainder between the index and the coordinate >*/
{
	/* declare remainder array */
	float* Rem = sf_floatalloc(ndim);
	/* declare looping index */
	int i ; 
    /* find coordinate for index, subtract from coordinate */
	float* Diff = conv_float_array_subtract( Coord, conv_index_to_coordinates( Ind, D, O, ndim), ndim);
	/* divide the difference by the increment */
	for ( i = 0 ; i < ndim ; i++){
		Rem [ i] = Diff[ i]/D[ i];
	}
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

int* conv_doughnut_wrap( int* Ind, int* N, int ndim)
	/*< checks to see if index is out of bounds.  If it is, performs doughnut wrapping >*/
{
	/* wrapped index */
	int* IndW = sf_intalloc(ndim);
	/* looping index */
	int i ;
	/* loop through dimensions */
	for ( i = 0 ; i < ndim ; i++ ){
		if ( Ind[ i] < 0 ){ 
			/* too low */
			IndW[ i] = Ind[ i] + N[ i];
		}else{
			if ( Ind[ i] >= N[ i]){
				/* too high */
				IndW[ i] = Ind[ i] - N[ i];
			}else{
				/* Goldilox */
				IndW[ i] = Ind[ i];
			}
		}
	}
	
	return IndW ;
}

float conv_interpolation_weights( int* NInd, float* Rem, int ndim)
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

float conv_array_doughnut_interpolator( int* Ind1, float* Rem, float* array, int* N, int ndim )
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
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* add that to the current index */
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* check to see if inbounds, if not, doughnut wrap */
		Ind3 = conv_doughnut_wrap( Ind2, N, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind2, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		interp += array[ indx] * conv_interpolation_weights( NInd, Rem, ndim);
	}
	return interp;
}

float conv_array_interpolator( int* Ind1, float* Rem, float* array, int* N, int ndim )
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
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* add that to the current index */
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind2, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		interp += array[ indx] * conv_interpolation_weights( NInd, Rem, ndim);
	}
	return interp;
}

float* conv_array_adj_interpolator( int* Ind1, float* Rem, float interp, float* arrayout, int* N, int ndim )
	/*< adjoint interpolates between elements of array >*/
{
	/* create N array for nodes */
	int* Nnode = conv_Nnode( ndim);
	/* Index array for node position */
	int* NInd = sf_intalloc( ndim);
	/* offset index array */
	int* Ind2 = sf_intalloc( ndim);
	/* how many nodes are we dealing with ? */
	long nnodes = (long)conv_int_exponent(2,ndim);
	/* looping index for nodes, using long for consistency with unwrapping programs */
	long nindx ;
	/* index for array position */
	long indx ;
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		NInd = conv_unwrap( nindx, Nnode, ndim);
		/* add that to the current index */
		Ind2 = conv_int_array_add( Ind1, NInd, ndim);
		/* unwrap to determine array position */
		indx = conv_wrap( Ind2, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		arrayout[ indx] += interp * conv_interpolation_weights( NInd, Rem, ndim);
	}
	return arrayout;
}

float* conv_translate(float* arrayin, float* X, int* N, float* D, float* O, int ndim)
	/*< returns a translated version of the array by vector X >*/	
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
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* check to see if we are in bounds */
		if ( conv_in_bounds( Ind2, N, ndim ) > 0 ) continue ; 
		/* interpolate, dont need += because looping through index */
		arrayout[ indx] = conv_array_interpolator( Ind2, TRem, arrayin, N, ndim );
	}
	/* return translated array */
	return arrayout;
}

float* conv_translate_wrap(float* arrayin, float* X, int* N, float* D, float* O, int ndim, bool adj)
	/*< returns a translated version of the array by vector X, with doughnut wrapping so an adjoint >*/	
{
	/* if adjoint, reverse translation */
	if (adj){
		X = conv_scale_float_array(X,-1.,ndim);
	}
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
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* interpolate */
	    arrayout[ indx] = conv_array_doughnut_interpolator( Ind2, TRem, arrayin, N, ndim ); 
	}
	/* return translated array */
	return arrayout;
}


float* conv_translate_adj_interp(float* arrayin, float* X, int* N, float* D, float* O, int ndim, bool adj)
	/*< returns a translated version of the array by vector X, with adjoint interpolation >*/	
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
		Ind1 = conv_unwrap( indx, N, ndim);
		/* Translate by X */
		Ind2 = conv_int_array_add( Ind1, TInd, ndim);
		/* check to see if we are in bounds */
		if ( conv_in_bounds(  Ind2, N, ndim ) > 0 ) continue ; 
		if(!adj){
			/* interpolate */
			 arrayout[ indx] = conv_array_interpolator( Ind2, TRem, arrayin, N, ndim ); }
		else{
			/* adjoint interpolation */ 
			arrayout = conv_array_adj_interpolator( Ind2, TRem, arrayin[ indx], arrayout, N, ndim ); }
	}
	/* return translated array */
	return arrayout;
}
