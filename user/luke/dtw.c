/* One Dimensional Dynamic Time Warping underlying programs, removing madagascar dependencies for arrays*/
/*^*/
#include <rsf.h>


void dtw_reverse_array(float* array, float* reverse, int n2, int n1){
	/*< reverses the second axis of an array >*/
	int i, j;
	for (i = 0 ; i < n2 ; i++){
		for ( j = 0 ; j < n1 ; j++){
			reverse[(n2-i-1)*n1 + j] = array[i*n1 + j] ;
		}
	}
	return;
}

void dtw_two_sided_smoothing_sum( float* ef, float* er, float* e, float* etilde, int n){
	/*< adds two arrays >*/
	int i;
	for ( i = 0 ; i < n ; i ++){
		if (ef[i] < 0 || er[i] < 0 || e[i] < 0) etilde[i] = -1;
		else etilde[i] = ef[i] + er[i] - e[i];
	}
	return;
}

float* _dtw_alignment_errors( float* match, float* ref, int n1, int maxshift , float ex){
	/*< finds pointwise alignment errors for shifts between min and max, returns mismatch >*/
	int j2, i, j;
	float* mismatch = malloc((n1*(2*maxshift+1))*sizeof(float));
	/* loop through samples in reference trace*/
	for ( i = 0 ; i < n1 ; i++){
		/* loop through permissible shifts */
		for ( j = 0 ; j < 2 * maxshift + 1 ; j++){
			j2 = i + ( j - maxshift );
			/* if in bounds write difference */
			if (j2 >= 0 && j2 < n1){ 
			    mismatch[ i*(2*maxshift+1) + j] = pow(fabsf(ref[ i] - match[ j2]),ex);
			}
			/* otherwise write a null value */
			else { mismatch [ i*(2*maxshift+1) + j] = -1.;}
		}
	}
	return mismatch;
}

void dtw_alignment_errors( float* match, float* ref, float* mismatch, int n1, int maxshift , float ex){
	/*< finds pointwise alignment errors for shifts between min and max >*/
	int j2, i, j;
	/* loop through samples in reference trace*/
	for ( i = 0 ; i < n1 ; i++){
		/* loop through permissible shifts */
		for ( j = 0 ; j < 2 * maxshift + 1 ; j++){
			j2 = i + ( j - maxshift );
			/* if in bounds write difference */
			if (j2 >= 0 && j2 < n1){ 
			    mismatch[ i*(2*maxshift+1) + j] = pow(fabsf(ref[ i] - match[ j2]),ex);
			}
			/* otherwise write a null value */
			else { mismatch [ i*(2*maxshift+1) + j] = -1.;}
		}
	}
	return;
}

float dtw_posmin(float a, float b){
	/*< finds the minimum non-null value if one exists >*/
	if( a >=0 && b >= 0){
		if ( a < b) return a;
		else return b;
	} else {
		if (a >= 0 ){
			return a;
		} else {
			return b;
		}
	}
}

float dtw_strain_accumulator(float* accumulate, float* mismatch, int i, int l, int s, int maxshift){
	/*< sums up bounded strains >*/
    /* initialize */
    float accum ;
	int j;
	if ( i-s > 0) accum	= accumulate[ (i-s)*(2*maxshift+1) + l];
	else accum = 0;
	for ( j = 0 ; j < s ; j++){
		if ( (i-j)*(2*maxshift+1) + l < 0 ) break ;
		if ( mismatch[(i-j)*(2*maxshift+1) + l] < 0) return mismatch[(i-j)*(2*maxshift+1) + l];
		accum += mismatch[(i-j)*(2*maxshift+1) + l];
	}
	return accum;
}

float dtw_acumin( float* accumulate, float* mismatch, int i, int l, int n1, int maxshift, int sb){
	/*< finds the minimum of 3 errors, checks for null values >*/
	/* read in the comparison values if in bounds */
	float first ;
	if (l-1 >=0){
	//    first = accumulate[ (i-1)*(2*maxshift+1) + l - 1];
        first = dtw_strain_accumulator(accumulate, mismatch, i, l-1, sb, maxshift);
	} else { 
		first = -1 ;
	}
	float second = accumulate[ (i-1)*(2*maxshift+1)+l];
	float third;
	if (l+1 < 2*maxshift+1 ){
	    third = dtw_strain_accumulator(accumulate, mismatch, i, l+1, sb, maxshift) ;
	} else { 
		third = -1 ;
	}	
	/* find the smallest non-null value */
	return dtw_posmin( first, dtw_posmin( second, third)) ;
}

void dtw_backtrack( float* accumulate, float* mismatch, int* shifts, int n1, int maxshift, float str){
	/*< backtrack to determine the shifts >*/
	int sb = (int)(1/str);
	/* determine initial shifts */
	int  u = maxshift;
	int du = 0;
	int l, s, i;
	float holder = 0;
	float test;
	for ( l = 0; l < 2*maxshift+1 ; l++){
		if ( accumulate [ (n1-1)*(2*maxshift+1) + l] < 0) continue;
		if ( accumulate [ (n1-1)*(2*maxshift+1) + l] < accumulate [ (n1-1)*(2*maxshift+1) + u]){
			u = l;
		}
	}
	/* write initial shift */
	shifts [ n1-1] = u - maxshift;	
	/* do initial shifts until we attain our strian bound */
	for ( s = 1 ; s < sb ; s++){
		u = shifts [n1-s];
		du = 0;
		/* holder = accumulate[ (n1-s-1)*(2*maxshift+1) + u + maxshift]; */
		holder = dtw_strain_accumulator(accumulate, mismatch, n1-s-1, u + maxshift, s, maxshift);
		if (u-1 + maxshift >= 0){
			test = dtw_strain_accumulator(accumulate, mismatch, n1-s-1, u-1 + maxshift, s, maxshift);
			if( test < holder ){
				holder = test;
				du = -1;
			}
		}
		if (u+1 + maxshift < 2*maxshift+1){
			test = dtw_strain_accumulator(accumulate, mismatch, n1-s-1, u+1 + maxshift, s, maxshift);
			if( test < holder ){
				holder = test;
				du = 1;
			}
		}
		shifts[ n1-s-1] = u+du;
	}
	/* loop backward through remaining time samples */
	for ( i = n1-sb-1 ; i >= 0 ; i--){
		u  = shifts[ i + 1 ];
		du = 0;
		holder = dtw_strain_accumulator(accumulate, mismatch, i, u + maxshift, sb, maxshift);
	    /* test above */
		if ( u-1 + maxshift >= 0 ){
			test = dtw_strain_accumulator(accumulate, mismatch, i, u-1 + maxshift, sb, maxshift);
			if (test < holder){
				holder = test;
				du = - 1;
			}
		}
		/* test below */
		if ( u+1 + maxshift < 2*maxshift+1 ){
			test = dtw_strain_accumulator(accumulate, mismatch, i, u+1 + maxshift, sb, maxshift);
			if (test < holder){
				du = 1;
			}
		}		
		/* write minimizing shift */
		shifts [ i] = u + du;
	}
	return;
}

void dtw_spread_over_nulls( float* accumulate, int n1, int maxshift){
	/*< spread the accumulation back over null values >*/
	int i, l;
	for (  i = 0; i < maxshift ; i++){
		for ( l = 0 ; l < maxshift - i; l++){
			/* beginning of signal */
			accumulate[i*(2*maxshift+1)+l] = 
				accumulate[(maxshift-l)*(2*maxshift+1)+l] ;
			/* end of signal */
			accumulate[(n1-i-1)*(2*maxshift+1)+(2*maxshift-l)] = 
				accumulate[(n1-1-maxshift+l)*(2*maxshift+1)+(2*maxshift-l)] ;
		}
	}
	return;
}

void dtw_accumulate_errors( float* mismatch, float* accumulate, int n1, int maxshift, float str){
	/*< accumulates errors over trajectories >*/
    /* strain bound */
	int sb = (int)(1/str);
	int s, l, i;
	/* initialize with first step */
	for ( s = 0 ; s < sb ; s++){
		for ( l = 0 ; l < 2*maxshift+1 ; l++){
			accumulate[s*(2*maxshift+1) + l] = dtw_strain_accumulator(accumulate, mismatch, s, l, s, maxshift);
		}
	}
	/* loop through rest of time indeces*/
	for ( i = sb ; i < n1 ; i++){
		/* loop through lags */
		for ( l = 0 ; l < 2*maxshift+1 ; l++){
			/* check to see if index makes sense, negative means null */
			if ( mismatch[ i*(2*maxshift+1) + l] < 0){
				accumulate[ i*(2*maxshift+1) + l ] = -1;
				continue;
			}
			/* Hale's equation 6 */
			float minval = dtw_acumin( accumulate, mismatch, i, l, n1, maxshift, sb);
			if (minval >= 0 ){
			    accumulate[ i*(2*maxshift+1) + l] = mismatch[ i*(2*maxshift+1) + l] + minval;
		    } else {
				/* write null val if all negative */ 
				accumulate[ i*(2*maxshift+1) + l] = minval ; 
			}
		}
	}
	return;
}

void dtw_apply_shifts( float* match, int* shifts, float* warped, int n){
	/*< apply known integer shifts to matching trace >*/
	int i;
	for ( i = 0 ; i < n ; i++){
		/* make sure shifts in bounds */
		if (i + shifts[ i] < 0){
			warped [ i] = match [ 0];
		} else { 
			if (i + shifts[ i] > n-1){
				warped[ i] = match [ n-1] ;
			} else{
				/* apply shifts */
				warped [ i] = match[ i + shifts[ i]]; 
			}
			
		}
	}
	return;
}

void dtw_apply_scaled_shifts( float* match, int* shifts, float* warped, float scale, int n){
	/*< apply scaled  integer shifts to matching trace >*/
	int i;
	for ( i = 0 ; i < n ; i++){
		/* read shift, scale */
		float shiftscl = scale*(float)shifts[ i];
		/* cast to integer */
		int shift = (int)shiftscl;
		/* make sure shifts in bounds */
		if (i + shift < 0){
			warped [ i] = match [ 0];
		} else { 
			if (i + shift > n-1){
				warped[ i] = match [ n-1] ;
			} else{
				/* apply shifts */
				warped [ i] = match[ i + shift]; 
			}
			
		}
	}
	return;
}

void dtw_copy(float* array, float val, int n)
/*< copies float vals into array >*/
{
	int i;
	for ( i = 0 ; i < n ; i++){
		array[i] = val;
	}
	return;
}


void dtw_acopy(float* array, float* vals, int n)
/*< copies vals into array, both size n >*/
{
	int i;
	for ( i = 0 ; i < n ; i++){
		array[i] = vals[i];
	}
	return;
}

void dtw_aadd(float* array, float* vals, int n)
/*< adds vals plus array, both size n >*/
{
	int i;
	for ( i = 0 ; i < n ; i++){
		array[i] += vals[i];
	}
	return;
}

void dtw_mul(float* array, float val, int n)
/*< scales array by float vals >*/
{
	int i;
	for ( i = 0 ; i < n ; i++){
		array[i] += val;
	}
	return;
}

void dtw_linear_comb(float* C, float alpha, float* A, float beta, float* B, int n)
	/*< C = alpha*A + beta*B >*/
{
	for ( int i = 0 ; i < n ; i++ ){
		C[ i] = alpha*A[ i] + beta*B[ i];
	}
	return;
}

void dtw_set_file_params(sf_file _file, int n1, float d1, float o1,
    const char* label1, const char* unit1,
    int n2, float d2, float o2,
	const char* label2, const char* unit2)
/*< set output parameters >*/	
{ 
		sf_putint   (_file,"n1",n1); 
		sf_putfloat (_file,"d1",d1);
		sf_putfloat (_file,"o1",o1);
		sf_putstring(_file,"label1", label1);
		sf_putstring(_file, "unit1",  unit1);
		sf_putint   (_file,"n2",n2);
		sf_putfloat (_file,"d2",d2);
		sf_putfloat (_file,"o2",o2);
		sf_putstring(_file,"label2", label2);
		sf_putstring(_file, "unit2",  unit2);
	return;
}

float* dtw_symmetric_accumulation(float* mismatch, int n1, int maxshift, float str)
	/*< performs symmetric smoothing accumulation given an input mismatch array, assumes we have already spread over nulls >*/
{
	float* accumulate_f = malloc(n1*(2*maxshift+1)*sizeof(float));
	/* accumulate forward errors */
	dtw_accumulate_errors( mismatch, accumulate_f, n1, maxshift, str);
	/* declare array for backward error accumulation */
	float* accumulate_b = malloc(sizeof(float)*n1*(2*maxshift+1));
	/* reverse mismatch */
	float* mismatch_r = malloc(sizeof(float)*n1*(2*maxshift+1));
	dtw_reverse_array(mismatch, mismatch_r, n1, (2*maxshift+1));
	/* accumulate backward errors */
	dtw_accumulate_errors( mismatch_r, accumulate_b, n1, maxshift, str);
	free (   mismatch_r);
	/* flip them */
	float* accumulate_br = malloc(sizeof(float)*n1*(2*maxshift+1));
	dtw_reverse_array(accumulate_b,accumulate_br, n1, (2*maxshift+1));
	free ( accumulate_b) ;
	/* sum the errors */
	float* accumulate = malloc(sizeof(float)*n1*(2*maxshift+1));
	dtw_two_sided_smoothing_sum(accumulate_f, accumulate_br, mismatch, accumulate,n1*(2*maxshift+1));
	free ( accumulate_f) ;
	free ( accumulate_br);
	
	return accumulate;
}

void dtw_symmetric_accumulation_v(float* accumulate, float* mismatch, int n1, int maxshift, float str)
	/*< performs symmetric smoothing accumulation given an input mismatch array, assumes we have already spread over nulls >*/
{
	float* accumulate_f = malloc(sizeof(float)*n1*(2*maxshift+1));
	/* accumulate forward errors */
	dtw_accumulate_errors( mismatch, accumulate_f, n1, maxshift, str);
	/* declare array for backward error accumulation */
	float* accumulate_b = malloc(sizeof(float)*n1*(2*maxshift+1));
	/* reverse mismatch */
	float* mismatch_r = malloc(sizeof(float)*n1*(2*maxshift+1));
	dtw_reverse_array(mismatch, mismatch_r, n1, (2*maxshift+1));
	/* accumulate backward errors */
	dtw_accumulate_errors( mismatch_r, accumulate_b, n1, maxshift, str);
	free (   mismatch_r);
	/* flip them */
	float* accumulate_br = malloc(sizeof(float)*n1*(2*maxshift+1));
	dtw_reverse_array(accumulate_b,accumulate_br, n1, (2*maxshift+1));
	free ( accumulate_b) ;
	/* zero out accumulation array */
	dtw_copy(accumulate, 0, n1*(2*maxshift+1));
	/* sum the errors */
	dtw_two_sided_smoothing_sum(accumulate_f, accumulate_br, mismatch, accumulate,n1*(2*maxshift+1));
	free ( accumulate_f) ;
	free ( accumulate_br);
	
	return;
}


void dtw_find_shifts(int* shifts, float* ref, float* match, 
    int n1, int maxshift, float str, float ex)
/*< integrated program to calculate shifts in one go, good for parallelization takes already allocated shifts with int(n1) as input>*/
{
	float* mismatch = malloc(sizeof(float)*n1*(2*maxshift+1));
	/* determine the best shifts using dynamic warping */
	dtw_alignment_errors( match, ref, mismatch, n1, maxshift, ex);
	/* spread values over nulls */
	dtw_spread_over_nulls( mismatch, n1, maxshift);
	/* declare array for forward error accumulation */
    float* accumulate = dtw_symmetric_accumulation( mismatch, n1, maxshift, str);
	/* backtrack to find integer shifts */
	dtw_backtrack( accumulate, mismatch, shifts, n1, maxshift, str);	
	/* remove unneded arrays */
	free ( accumulate );
	free ( mismatch );
    return;
}

float* dtw_write_gather_nullvals(float* gather, int n1, int n2, float nullval)
	/*< writes null values to gatherarray >*/
{
	float* ngather = malloc(sizeof(float)*n1*n2);
	dtw_acopy(ngather, gather, n1*n2);
	int i1, i2;
	/* loop through traces */
	for ( i2 = 0; i2 < n2 ; i2++ ){
		/* loop down */
		for ( i1 = 0 ; i1 < n1 ; i1++){
			if (ngather[ i1 + i2*n1] == 0.){
				/* write null value */
				ngather [ i1 + i2*n1] = nullval;
			} else { break ;}
		}
		/* loop up */
		for (i1 = n1-1 ; i1 >= 0 ; i1--){
			if (ngather[ i1 + i2*n1] == 0.){
				/* write null value */
				ngather [ i1 + i2*n1] = nullval;
			} else { break ;}
		}
		
	}
	return ngather;
}

float* dtw_copy_nulls(float* clean, float* nulls, int n, float nullval)
	/*< puts a nullval in the clean array with the same index if one exists in the nulls array >*/
{
	float* nclean = malloc(sizeof(float)*n);
	int i;
	for (i = 0; i < n ; i++){
		if ( nulls[ i] == nullval){
			/* write null to clean */
			nclean[ i] = nullval;
		} else{
			nclean[ i] = clean [ i];
		}
	}
	return nclean;
}

float* dtw_subtract_minimum(float* array, int n)
	/*<find minimum value of float array and subract it >*/
{
	int i;
	float minval = SF_FLOAT;
	float* arrout = malloc(sizeof(float)*n);
	/* find minimum */
	for ( i = 0 ; i < n ; i++ ){
		if (array[ i] < minval){
			minval = array[ i];
		}
	}
	/* subtract minimum */
	for ( i = 0 ; i < n ; i++){
		arrout[ i] = array[ i] - minval;
	}
	return arrout;
}

void dtw_subtract_minimum_v(float* arrout, float* array, int n)
	/*<find minimum value of float array and subract it >*/
{
	int i;
	float minval = SF_FLOAT;
	/* find minimum */
	for ( i = 0 ; i < n ; i++ ){
		if (array[ i] < minval){
			minval = array[ i];
		}
	}
	/* subtract minimum */
	for ( i = 0 ; i < n ; i++){
		arrout[ i] = array[ i] - minval;
	}
	return ;
}

void dtw_norm_stack(float* gather, float* stack, int n1, int n2, float nullval)
/*< normalized stacking >*/
{
	/* clear out stack array */
    dtw_copy( stack, 0., n1);
	/* loop indexes */
	int i1, i2;
	/* fold array */
	int*    fold = malloc(sizeof(int)*n1);
	/* loop through */
	for ( i2 = 0 ; i2 < n2 ; i2++ ){
		for ( i1 = 0 ; i1 < n1 ; i1++ ){
			if ( gather[ i2*n1 + i1 ] == nullval) continue ; 
			stack[ i1] += gather [ i2*n1 +i1 ];
			fold [ i1] += 1;
		}
	}
	/* normalize by fold */
	for ( i1 = 0 ; i1 < n1 ; i1++ ){
		if ( fold[ i1] == 0 ) { stack [ i1] = 0.; }
		else { stack [ i1 ] = stack [ i1 ] / ( (int) fold[ i1 ] ); }
	}
	free(fold);
	return;
}

void dtw_get_column( float* array, float* column, int i, int n )
/*< grab ith column from an array with height n>*/	
{
	
	int j;
	for ( j =  0 ; j < n ; j++ ){
		column [ j ] = array [ j + i*n ];
	}
	return;
}

float* _dtw_get_column( float* array, int i, int n )
/*< grab ith column from an array with height, outputs the column n>*/	
{
	float* column = malloc(sizeof(float)*n);
	int j;
	for ( j =  0 ; j < n ; j++ ){
		column [ j ] = array [ j + i*n ];
	}
	return column;
}

float* dtw_int_to_float( int* intarray, int n)
/*< write an int array to a float array >*/
{
	float* floatarray = malloc(sizeof(float)*n);
	int i;
	for ( i = 0 ; i < n ; i++){
		floatarray[ i] = (float)intarray[ i];
	}
	return floatarray;
}

void dtw_int_to_float_v( float* floatarray, int* intarray, int n)
	/*< write an int array to a float array >*/
{
	for (int i = 0 ; i < n ; i++){
		floatarray[i] = (float)intarray[i];
	}
	return;
}

void dtw_float_to_int_v( int* intarray, float* floatarray, int n)
	/*< write a float array to an int array >*/
{
	for (int i = 0 ; i < n ; i++ ){
		intarray[i] = (int)floatarray[i];
	}
	return;
}

void dtw_put_column( float* array, float* column, int i, int n )
/*< put ith column into an array with height n>*/	
{
	
	int j;
	for ( j =  0 ; j < n ; j++ ){
		array [ j + i*n ] = column [ j ] ;
	}
	return;
}

long dtw_size(int* N, int ndim)
/*< findes size of an array by multiplying out dimensions >*/
{
	long sz = 1;
	int i;
	for (i = 0 ; i < ndim ; i++){
		sz *= (long)N[ i];
	}
	return sz;
}


int* dtw_unwrap(long indx, int* N, int ndim)
	/*< unwraps index into coordinates >*/
{
    int* Ind = malloc(sizeof(int)*ndim);
	int j, n;
	long r = indx;
	for (j = ndim-1 ; j >= 0 ; j--){
        n = dtw_size(N,j);
		Ind [ j] = (int)((double)r / (double)n);
		r -= (long) Ind[j]*n;
	}
	return Ind;
}

long dtw_wrap(int* Ind, int* N, int ndim)
/*< wraps input coordinates back into index >*/
{
	long indx = 0;
	int j, n; 
	for ( j = 0; j < ndim ; j++){
        n = dtw_size(N,j);
		indx += n * Ind[j];
	}
	return indx;
}


int* dtw_swaperoo(int* In, int t1, int t2, int ndim)
	/*< swaps the t1 and t2 of In >*/
{
	int* Out = malloc(sizeof(int)*ndim);
	int i ;
	/* initialize */
	for (i = 0 ; i < ndim ; i++){
		Out [ i] = In [ i];
	}
	/* CHANGE PLACES! */
	int holdr1 = Out [ t1 - 1];
	int holdr2 = Out [ t2 - 1];
	Out [ t2 - 1 ] = holdr1;
	Out [ t1 - 1 ] = holdr2;
	return Out;
}


void dtw_swaperoo_v(int* Out, int* In, int t1, int t2, int ndim)
	/*< swaps the t1 and t2 of In >*/
{
	int i ;
	/* initialize */
	for (i = 0 ; i < ndim ; i++){
		Out [ i] = In [ i];
	}
	/* CHANGE PLACES! */
	int holdr1 = Out [ t1 - 1];
	int holdr2 = Out [ t2 - 1];
	Out [ t2 - 1 ] = holdr1;
	Out [ t1 - 1 ] = holdr2;
	return ;
}

float* dtw_swaperoof(float* In, int t1, int t2, int ndim)
	/*< swaps the t1 and t2 of In for floats >*/
{
	float* Out = malloc(sizeof(float)*ndim);
	int i ;
	/* initialize */
	for (i = 0 ; i < ndim ; i++){
		Out [ i] = In [ i];
	}
	/* CHANGE PLACES! */
	float holdr1 = Out [ t1 - 1];
	float holdr2 = Out [ t2 - 1];
	Out [ t2 - 1 ] = holdr1;
	Out [ t1 - 1 ] = holdr2;
	return Out;
}

void dtw_swaperoof_v(float* Out, float* In, int t1, int t2, int ndim)
	/*< swaps the t1 and t2 of In for floats >*/
{
	int i ;
	/* initialize */
	for (i = 0 ; i < ndim ; i++){
		Out [ i] = In [ i];
	}
	/* CHANGE PLACES! */
	float holdr1 = Out [ t1 - 1];
	float holdr2 = Out [ t2 - 1];
	Out [ t2 - 1 ] = holdr1;
	Out [ t1 - 1 ] = holdr2;
	return ;
}

float* dtw_transp( float* arr, int t1, int t2, int* Nin, int ndim)
	/*< follows madagascar convention for transp plane=t1t2 >*/
{
	/* allocate transpose array */
	float* arrT = malloc(sizeof(float)*dtw_size(Nin,ndim));
	/* determine N of output array */
	int* Nout = dtw_swaperoo( Nin, t1, t2, ndim);
	/* loop through array elements and transpose */
	long i;
	for (i = 0 ; i < dtw_size(Nin,ndim); i++){
		arrT[ dtw_wrap(dtw_swaperoo(dtw_unwrap(i, Nin, ndim),t1,t2,ndim),Nout,ndim)] = arr [ i];
	}
	return arrT;
	free (Nout);
}

void dtw_transp_v(float* arrT, float* arr, int t1, int t2, int* Nin, int ndim)
	/*< follows madagascar convention for transp plane=t1t2 >*/
{
	/* determine N of output array */
	int* Nout = dtw_swaperoo( Nin, t1, t2, ndim);
	/* index */
	int* Ind;
	/* loop through array elements and transpose */
	long i;
	for (i = 0 ; i < dtw_size(Nin,ndim); i++){
		/* find index */
		Ind = dtw_unwrap(i, Nin, ndim);
		/* swap */
		dtw_swaperoo_v(Ind,Ind,t1,t2,ndim);
		arrT[ dtw_wrap(Ind,Nout,ndim)] = arr [ i];
		/* remove dynamic array */
		free ( Ind);
	}
	/* free dynamic array */
	free (Nout);
	return;
}
