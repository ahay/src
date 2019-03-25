/* Path Integral Minimization Tools */
/*^*/
#include <rsf.h>
#include <math.h>

float path_abs( float a)
	/*< returns absolute value of a >*/
{
	if ( a < 0 ) a *= -1;
	return a;
}

float path_max( float a, float b)
	/*< returns max of a, b >*/
{
	if ( a > b ) return a;
	else return b;
}

float path_min( float a, float b)
	/*< returns min of a, b >*/
{
	if ( a > b ) return b;
	else return a;
}

long path_size(int* N, int ndim)
/*< findes size of an array by multiplying out dimensions >*/
{
	long sz = 1;
	int i;
	for (i = 0 ; i < ndim ; i++){
		sz *= (long)N[ i];
	}
	return sz;
}


void path_unwrap(int* Ind, long indx, int* N, int ndim)
	/*< unwraps index into coordinates >*/
{
	int j, n;
	long r = indx;
	for (j = ndim-1 ; j >= 0 ; j--){
        n = path_size(N,j);
		Ind [ j] = (int)((double)r / (double)n);
		r -= (long) Ind[j]*n;
	}
	return;
}


long path_wrap(int* Ind, int* N, int ndim)
/*< wraps input coordinates back into index >*/
{
	long indx = 0;
	int j, n; 
	for ( j = 0; j < ndim ; j++){
        n = path_size(N,j);
		indx += n * Ind[j];
	}
	return indx;
}

void path_swaperoof(float* Out, float* In, int t1, int t2, int ndim)
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


void path_swaperoo(int* Out, int* In, int t1, int t2, int ndim)
	/*< swaps the t1 and t2 of Int arrays for indicies >*/
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


void path_transp(float* arrT, float* arr, int t1, int t2, int* Nin, int ndim)
	/*< follows madagascar convention for transp plane=t1t2 >*/
{
	/* determine N of output array */
	int* Nout = sf_intalloc(ndim);
	path_swaperoo(Nout, Nin, t1, t2, ndim);
	/* index */
	int* Ind = sf_intalloc(ndim);
	/* loop through array elements and transpose */
	long i;
	for (i = 0 ; i < path_size(Nin,ndim); i++){
		/* find index */
		path_unwrap(Ind, i, Nin, ndim);
		/* swap */
		path_swaperoo(Ind,Ind,t1,t2,ndim);
		arrT[ path_wrap(Ind,Nout,ndim)] = arr [ i];
		/* remove dynamic array */
	}
	/* free dynamic array */
	free (Nout);
	free (Ind);
	return;
}


void path_combine(float* C, float* A, float a, float* B, float b, long n)
	/*< C = aA + bB >*/
{
	/* looping index */
	long i ;
	for ( i = 0 ; i < n ; i++ ){
		C [ i] = a*A[i] + b*B[i] ;
	}
	return;
}

void path_scale(float* C, float* A, float a, long n)
	/*< C = aA >*/
{
	/* looping index */
	long i ;
	for ( i = 0 ; i < n ; i++ ){
		C [ i] = A [ i] * a ; 
	}
	return;
}

void path_copy_i( int* C, int* A, long n)
	/*< copy A into C >*/
{
	/* looping index */
	long i;
	for ( i = 0 ; i < n ; i++){
		C[ i] = A[ i];
	}
	return;
}

void path_subtract( float* C, float* A, float* B, long n)
	/*< C = A - B >*/
{
	path_combine(C, A, 1., B, -1., n);
	return;
}


void path_add_i( int* C, int* A, int* B, long n)
	/*< arr int arrays A+B=C >*/
{
	/* looping index */
	long i;
	for (i = 0 ; i < n ; i++){
		C[ i] = A [i]+B[i];
	}
	return;
}
	

float path_dot( float* A, float* B, long n)
	/*< cannonical dot product of A dot B >*/
{
	/* looping index */
	long i;
	/* initialize dot */
	double dot = 0;
	for ( i = 0 ; i < n ; i++){
		dot += A[ i] * B[ i];
	}
	return (float) dot;
}

float path_norm ( float* A, long n)
	/*< calculates L2 norm of A >*/
{
	return sqrtf(path_dot(A,A,n));
}

float path_distance(float* A, float* B, long n)
	/*< computes l2 distance between arrays A, B >*/
{
	/* array for holding the difference */	
	float* dAB = sf_floatalloc(n);
	/* compute difference */
	path_subtract(dAB,A,B,n);
	/* compute distance */
	float dist = path_norm(dAB,n);
	/* clear intermediate array */
	free(dAB);
	/* return distance */
	return dist;
}

void path_normalize ( float* N, float* A, long n)
	/*< normalizez vector A of length n >*/
{
	/* magnitude */
	float mag = path_norm( A, n);
	if ( mag > 0 ){
		/* avoid division by zero */
		path_scale( N, A, 1./mag, n);
	} else {
		/* write zeroes */
		path_scale( N, A, 0, n);
	}
	return;
}

void path_zero( float* A, long n)
	/*< zeros out a n array >*/
{
	path_scale( A, A, 0, n);
	return;
}

void path_zero_i( int* A, long n)
	/*< zeros out a n array for integers >*/
{
    long i ;
	for (i = 0 ; i <n ; i++ ){
		A[ i] = 0;
	}
	return;
}

void path_project(float* C, float* A, float* B, long n)
	/*< C =   A projected onto range of B >*/
{
	/* normalize B */
	float* B_norm = sf_floatalloc(n);
	path_normalize( B_norm, B, n);
	/* looping index */
	long i, j;
	/* initialize C */
	path_zero(C,n);
	/* project */
	for ( i = 0 ; i < n ; i++ ){
		for ( j = 0 ; j < n ; j++){
			C[ i] += A [ j] * B_norm[ j];
		}
		C [ i] *= B_norm[ i];
	}
	/* free memory */
	free (B_norm);
	return;
	
}

void path_orthagonalize(float* C, float* A, float* B, long n)
	/*< C =   component of A perpendicular to B >*/
{
	/* normalize B */
	float* B_norm = sf_floatalloc(n);
	path_normalize( B_norm, B, n);
	/* initialize C */
	path_zero(C,n);
	/* project */
    path_project(C, A, B, n);
	/* get perpindicular component */
	path_combine(C,A,1,C,-1,n);
	/* free memory */
	free (B_norm);
	return;
	
}

void path_get_column( float* array, float* column, int i, int n )
/*< grab ith column from an array with height n>*/	
{
	
	int j;
	for ( j =  0 ; j < n ; j++ ){
		column [ j ] = array [ j + i*n ];
	}
	return;
}

void path_put_column( float* array, float* column, int i, int n )
/*< put ith column into an array with height n>*/	
{
	
	int j;
	for ( j =  0 ; j < n ; j++ ){
		array [ j + i*n ] = column [ j ] ;
	}
	return;
}

void path_make_max( float* Max, int* N, float* D, float* O, int ndim)
	/*< make min and max arrays >*/
{
	/* looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++ ){
		Max[ i] = ((float)N[i]-1) * D[i] + O[i];
	}
	return;
}

void path_resample_arrays( float* C ,float* A, float* B, int nnew_knots, int ndim)
	/*< resample between A & B into C with N extra points in between >*/
{
	/* total knots */
	int nknots = nnew_knots + 2;
	/* fraction of distance along path */
	float frac;
	/* looping index */
	int i;
	/* zero out C */
	path_zero( C,  (long)ndim*(nnew_knots+2));
	/* holder coordinate array */
	float* Coord = sf_floatalloc(ndim);
	/* loop through knots */
	for ( i = 0 ; i < nknots ; i++){
		/* how far alont are we */
		frac = (float)i / (float) (nknots - 1);
		/* interpolate */
		path_combine( Coord, A, ( 1.-frac ), B, frac, ndim);
		/* put in output array */
		path_put_column( C, Coord, i, ndim );
	}
	/* free coordinate array */
	free (Coord);
	return;
}

void path_one_dim_smooth(float* Smooth, float* In, int order, long n)
	/*< one dimensional smoothing >*/
{
	/* looping index */
	long i, counter, io, this;
	float holder;
	for ( i = 0 ; i < n ; i++){
		counter = 0;
		holder = 0;
		for ( io = -1*order; io <= order ; io++){
			this = io + i;
			if ( this < 0 || this > n ) continue;
			holder += In[this];
			counter +=1;
		}
		if (counter > 0){
			Smooth[ i] = holder / (float)counter;
		}else{
			Smooth[i] = 0.;
		}
	}
	return;
}


void path_deriv(float* dA, float* A, int n1, int n2, float d1, int dorder, int slen, int nsmooth, bool scale)
	/*< calculate derivative of array A in first direction >*/
{
	/* initiate derivative */
    sf_deriv_init(n1, dorder, 0.);
	/* allocate data */
	float* dat = sf_floatalloc(n1);
	/* and its derivative */
	float* der = sf_floatalloc(n1);
	/* looping idexes */
	int i1, i2;
    for (i2=0 ; i2 < n2 ; i2++ ) {
		/* get data vector */
		path_get_column(A,dat,i2, n1);
		sf_deriv(dat,der);
		/* scale */
		if (scale){
			for (i1=0; i1 < n1; i1++) {
				der[i1] /= d1;
			}
		}
		/* smoothing */
		for ( i1 = 0 ; i1 < nsmooth ; i1++){
			path_one_dim_smooth(der, der, slen, n1);
		}
		/* put derivative in array */
		path_put_column(dA,der,i2, n1);
	}

	/* free arrays */
	free (dat );
	free (der);
	return;
}

void path_first_guess1(float* Ro, int knots, int* N, float* D, float* O, int ndim)
	/*< generate our first guess R >*/
{
	/* find array maximum */
	float* M  = sf_floatalloc(ndim);
	path_make_max(M , N , D, O, ndim);
	/* create fist point */
	float* Xo = sf_floatalloc(ndim);
	path_zero(Xo,ndim);
	Xo [ ndim-1] = (M [ ndim-1] + O[ ndim-1])/2;
	
	/* and last */
	float* Xf = sf_floatalloc(ndim);
	int i;
	for ( i = 0 ; i < ndim ; i++ ){
		Xf[ i] = M [ i];
	}
    Xf [ ndim-1] = (M [ ndim-1] + O[ ndim-1])/2;
	/* interpolate into path */
	path_resample_arrays( Ro , Xo, Xf, knots-2, ndim);
	/* free arrays */
	free (M);
	free( Xo);
	free( Xf);
	return;
}

void path_path_derivative(float* dR, float* R,int knots, int ndim)
/*< calculeate dRi = Ri - R(i-1) >*/
{
	/* zero out derivative */
	path_zero(dR, (knots+2)*ndim);
	/* looping index */
	int i, j;
	for (j = 1 ; j < knots ; j++){
		for ( i = 0 ; i < ndim ; i++ ){
			dR[ j*ndim + i ] = R [ j*ndim + i ] - R [ (j-1)*ndim + i ];
		}
	}
	return;
}

void path_index_to_coordinates(float* Coord, int *Ind, float *D, float *O, int ndim)
	/*< transforms an index location to a coordinate >*/
{
	/* declare looping index */
	int i ;
	for (i = 0 ; i < ndim ; i++){
		Coord[ i] = ((float) Ind[ i]) * D[ i] + O[ i];
	}
	return ;
}



void path_coordinates_to_index(int* Ind, float *Coord, float *D, float *O, int ndim)
	/*< converts coordinate array to nearest index rounding down >*/
{
	/* looping index */
	int i ;
	/* loop through dimensions */
	for ( i = 0 ; i < ndim ; i++){
		Ind[ i] = (int)((Coord[ i] - O[ i])/D[i]);
	}
	return ;
}

void path_index_coords_remainder(float* Rem, int *Ind, float *Coord, float *D, float *O, int ndim)
	/*< finds the interpolation proportion remainder between the index and the coordinate >*/
{
	/* declare converted coord array */
    float* CoordC = sf_floatalloc(ndim);
	/* get converted coordinates */
	path_index_to_coordinates(CoordC, Ind, D, O, ndim);
	/* declare looping index */
	int i ; 
    /* find coordinate for index, subtract from coordinate */
	float* Diff = sf_floatalloc(ndim);
	path_combine(Diff, Coord,1., CoordC,-1, ndim);
	/* divide the difference by the increment */
	for ( i = 0 ; i < ndim ; i++){
		Rem [ i] = Diff[ i]/D[ i];
	}
	free (Diff);
	free (CoordC);
}



void path_Nnode(int* Nnode, int ndim)
	/*< returns a N array for the nodes where each n is 2 >*/
{
	/* declare looping index */
	int i;
	for ( i = 0 ; i < ndim ; i++){
		Nnode[ i] = 2;
	}
	return;
}

float path_interpolation_weights( int *NInd, float *Rem, int ndim)
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

float path_array_interpolator( int *Ind1, float *Rem, float *array, int *N, int ndim )
	/*< interpolates between elements of array >*/
{
	/* value to be interpolated */
	float interp = 0;
	/* create N array for nodes */
	int* Nnode = sf_intalloc(ndim);
	path_Nnode(Nnode, ndim);
	/* Index array for node position */
	int* NInd = sf_intalloc( ndim);
	/* offset index array */
	int* Ind2 = sf_intalloc( ndim);
	/* how many nodes are we dealing with ? */
	long nnodes = path_size( Nnode, ndim);
	/* looping index for nodes, using long for consistency with unwrapping programs */
	long nindx ;
	/* index for array position */
	long indx ;
	/* loop through node points for interpolation */
	for ( nindx = 0 ; nindx < nnodes ; nindx++ ){
		/* get where we are */
		path_unwrap(NInd, nindx, Nnode, ndim);
		/* add that to the current index */
		path_add_i(Ind2, Ind1, NInd, ndim);
		/* unwrap to determine array position */
		indx = path_wrap( Ind2, N, ndim);
		/* read that value from array, add weighted value to interpolation */
		interp += array[ indx] * path_interpolation_weights( NInd, Rem, ndim);
	}
	/* free unneeded arrays */
	free (Nnode);
	free (NInd);
	free (Ind2);
	return interp;
}

int path_in_bounds_interp( int *Ind, int *N, int ndim )
	/*< determines if an index is in bounds, return 0 if in bounds, 1 if out of bounds >*/
{
	/* looping index */
	int i;
	/* loop through to see if we go out of bounds */
	for ( i = 0 ; i < ndim ; i++ ){
		if ( Ind[ i] < 0 ) return 1;
		if ( Ind[ i] > N[ i] - 1 ) return 1;
		if ( Ind[ i] == N[ i] - 1) return 2+i;
	}
	return 0;
}

float path_evaluate_array(float* X, float* Array, int* N, float* D, float* O, int ndim)
	/*< evaluate Array at X >*/
{
	/* index */
	int* Ind = sf_intalloc(ndim);
	float* Rem = sf_floatalloc(ndim);
	/* convert X to index */
	path_coordinates_to_index(Ind, X, D, O, ndim);
	/* and get the remainder */
	path_index_coords_remainder(Rem, Ind, X, D, O, ndim);
	/* interpolate */
	float interp = 0 ;
	int inbounds = path_in_bounds_interp(Ind, N, ndim) ;
	if ( inbounds == 0  ){ 
		interp = path_array_interpolator( Ind, Rem, Array, N, ndim );
	}else{
		if (inbounds == 1){
			sf_error("Interpolator going out of bounds! %i",path_in_bounds_interp(Ind, N, ndim));
		} else {
			inbounds -= 2;
			Ind [ inbounds] -=1;
			Rem [ inbounds]  = 1;
			interp = path_array_interpolator( Ind, Rem, Array, N, ndim );
		}
	}
	/* free arrays */
	free (Ind);
	free (Rem);
	return interp ;
}

void path_grad_s1(float* G,float* dS, int* N, float* D, float* O, float* R, float* dR, int knots, int ndim)
	/*< calculates the gradient of S for 1d case >*/
{
	/* zero out gradient */
	path_scale(G,G,0,knots*(ndim+1));
	/* knot index */
	int ik;

	/* dimension index */
	int id;
	/* distance element*/
	float dscale;
	/* coordinate array */
	float* X = sf_floatalloc(ndim);
	/* and the DRs */
	float* dRo = sf_floatalloc(ndim);
	float* dRf = sf_floatalloc(ndim);
	/* make dummy dimensions */
	int* N1 = sf_intalloc(ndim+1);
	float* D1 = sf_floatalloc(ndim+1);
	float* O1 = sf_floatalloc(ndim+1);

	path_scale(D1, D, 1, ndim);
	D1[ndim] = 1;
	path_scale(O1, O, 1, ndim);
	O1[ndim] = 0;
	path_copy_i(N1, N, ndim);
	N1[ndim] = ndim;

	float* X1 = sf_floatalloc(ndim+1);
	/* loop through knots */
	for ( ik = 0 ; ik < knots ; ik++){
		/* get where we are */
		path_get_column(R, X, ik, ndim);
		/* and the dRs */
		path_get_column(dR, dRo, ik, ndim);
		path_get_column(dR, dRf, ik+1, ndim);
		/* distance */
		dscale = (path_norm ( dRo, ndim) + path_norm ( dRf, ndim));
		/* loop through dimensions */
		for ( id = 0 ; id < ndim ; id++){
			/* clear X1 */
			path_scale(X1,X1,0,ndim+1);
			/* write dummy amounts */
			path_scale(X1, X, 1, ndim);
			X1 [ndim ] = (float)id*D1[ndim];
		    /* ealuate dS at this point, we are only allowing this to perturb  in first dimension*/
		    G [ ik*ndim + id ] = path_evaluate_array(X1, dS, N1, D1, O1 , ndim+1) * dscale ; 
		}
	}
	/* free arrays */
	free(X);
	free(X1);
	free (dRo);
	free (dRf);
	free( N1);
	free( O1);
	free (D1);
	
	return;
}

void path_gradient(float* Sfunc,float* dSfunc, int* N, float* D, float* O, int dorder, int slen, int nsmooth, int ndim, bool scale)
	/*< get gradient >*/
{
	long panelsize = path_size(N,ndim);
	float* St = sf_floatalloc(panelsize);
	float* dSt = sf_floatalloc(panelsize);
	float* dS = sf_floatalloc(panelsize);
	/* get derivative */
	int* Nt = sf_intalloc(ndim);
	float* Dt = sf_floatalloc(ndim);

	int dim, t1, t2;
/*	
	int t11, t22;
	*/
	for ( dim = 0 ; dim < ndim ; dim++){
		if ( dim > 0){
		   t1 = 1;
		   t2 = dim+1;
		   path_transp(St, Sfunc, t1, t2, N, ndim);
		   path_swaperoo(Nt, N, t1, t2, ndim);
		   path_swaperoof(Dt, D, t1, t2, ndim);
	   }else{
		   /* copy array */
		   path_scale(St, Sfunc, 1, panelsize);
		   path_copy_i(Nt,N,ndim);
		   path_scale(Dt, D, 1, ndim);
	   }
/*	   if (ndim > 2){
		   t11 = 2;
		   t22 = ndim ;
		   path_transp(St, St, t11, t11, Nt, ndim);
		   path_swaperoo(Nt, Nt, t11, t22, ndim);
		   path_swaperoof(Dt, Dt, t11, t22, ndim);
	   }
	   */
	   	path_deriv(dSt,St,Nt[0],panelsize/Nt[0],Dt[0],dorder, slen, nsmooth, scale);
		/* and go back */
/*	   if (ndim > 2){
		   t11 = 2;
		   t22 = ndim;
		   path_transp(dSt, dSt, t11, t22, Nt, ndim);
		   path_swaperoo(Nt, Nt, t11, t22, ndim);
	   }
		*/
	   if ( dim > 0){
	       path_transp(dS, dSt, t1, t2, Nt, ndim);
	   } else {
	   	   path_scale(dS, St, 1, panelsize);
	   }
	   /* write in the gradient array */
	   path_put_column(dSfunc, dS, dim, panelsize);
	}
	free (Nt);
	free (Dt);
	free(dSt);
	free (St);
	free (dS);
}


void path_create_tau_plus_minus(float* Tau_plus, float* Tau_minus, float* R, int nknots, int ndim)
	/*<create the Tau approximation by averaging dR around point>*/
{
	/* looping index */
	int i;
	/* local arrays */
	float* Tau_mL = sf_floatalloc(ndim);
	float* Tau_pL = sf_floatalloc(ndim);
	float* R0  = sf_floatalloc(ndim);
	float* R1  = sf_floatalloc(ndim);
	float* R2  = sf_floatalloc(ndim);
	float d1;
	float d2;
	float* D2 = sf_floatalloc(ndim);
	bool last = false;
	/* zero out tau */
	path_scale( Tau_plus, Tau_plus, 0, nknots*ndim);
	path_scale( Tau_minus, Tau_minus, 0, nknots*ndim);
	for ( i = 0 ; i < nknots ; i++ ){
		/* get comuns */
		/* centered */
		path_get_column(R, R1, i, ndim);
		/* next */
		if ( i < nknots - 1) {
			path_get_column(R, R2, i+1, ndim);
		} else {
			/* avoid pulling a null */
			path_scale(R2, R2, 0, ndim);
			last = true;
		}
		/* subtract */

		/* previous */
		if ( i > 0) {
			path_get_column(R, R0, i-1, ndim);
		} else {
			/* extrapolate */
			/* really R3 */
			path_get_column(R, R0, i+2, ndim);
			/* difference */
			path_combine(R0, R2, -1, R0, 1, ndim);
			/* get the other tau for 2nd deriv */
			path_combine(Tau_pL, R1, -1.0, R2, 1.0, ndim);
			/* distances */
			d2 = path_norm(R0,ndim);
			d1 = path_norm(Tau_pL,ndim);
			/* second derivative */
			path_combine(D2, R0, 2*d2/(d2+d1), Tau_pL, -2*d1/(d2+d1), ndim);
			/* apply */
			path_combine(R0, R0, 1, D2, (d2*d2+d1*d1)/(d1*d1), ndim);
			/* project backward */
			path_combine(R0, R0, -(d1+d2)/d1, R2, 1, ndim);
		}		
		/* worry about the last index */
		if ( last ){
			/* extrapolate */
			/* really R-1 */
			path_get_column(R, R2, i-2, ndim);
			/* difference */
			path_combine(R2, R2, -1, R0, 1, ndim);
			/* get the other tau for 2nd deriv */
			path_combine(Tau_mL, R0, -1.0, R1, 1.0, ndim);
			/* distances */
			d2 = path_norm(R2,ndim);
			d1 = path_norm(Tau_mL,ndim);
			/* second derivative */
			path_combine(D2, R2, 2*d2/(d2+d1), Tau_mL, -2*d1/(d2+d1), ndim);
			/* apply */
			path_combine(R2, R2, 1, D2, (d2*d2+d1*d1)/(d1*d1), ndim);
			/* project forward */
			path_combine(R2, R2, (d1+d2)/d1, R1, 1, ndim);			
		}
		/* subtract */
		path_combine(Tau_pL, R1, -1.0, R2, 1.0, ndim);
		path_combine(Tau_mL, R0, -1.0, R1, 1.0, ndim);

		/* put column back */
		path_put_column(Tau_plus,Tau_pL,i,ndim);
		path_put_column( Tau_minus, Tau_mL, i, ndim);
	}
	/* free arrays */
    	free (Tau_mL);
		free (Tau_pL);
		free ( R0);
		free ( R1);
		free ( R2);
		free (D2);
	return;
}

void path_evaluate_potental(float* V, float* R, int nknots, float* Potential, int* N, float* D, float* O, int ndim)
	/*< evaluate the potentail at the knot points >*/
{
	/* looping index */
	int i ;
	/* local array */
	float* X  = sf_floatalloc(ndim);
	/* clear out potential */
	path_scale(V,V,0., nknots);
	/* loop through knots */
	for ( i = 0 ; i < nknots ; i++ ){
		/* get column */
		path_get_column(R, X, i, ndim);
		/* evaluate potential */
		V[i] = path_evaluate_array(X, Potential, N, D, O , ndim);
	}
	/* free intermediate array */
	free (X);
	return;
}

void path_evaluate_gradient(float* G, float* R, int nknots, float* dPotential, int* N, float* D, float* O, int ndim)
	/*< evaluate the gradient of our objective function at knots defininng R >*/
{
	/* zero out gradient */
	path_scale(G,G,0,nknots*ndim);
	/* knot index */
	int ik;
	/* dimension index */
	int id;
	/* coordinate array */
	float* X = sf_floatalloc(ndim);
	/* make dummy dimensions */
	int* N1 = sf_intalloc(ndim+1);
	float* D1 = sf_floatalloc(ndim+1);
	float* O1 = sf_floatalloc(ndim+1);
	
	/* create dummy dimension arrays */
	path_scale(D1, D, 1, ndim);
	D1[ndim] = 1;
	path_scale(O1, O, 1, ndim);
	O1[ndim] = 0;
	path_copy_i(N1, N, ndim);
	N1[ndim] = ndim;

	float* X1 = sf_floatalloc(ndim+1);
	/* loop through knots */
	for ( ik = 0 ; ik < nknots ; ik++){
		/* get where we are */
		path_get_column(R, X, ik, ndim);
		/* loop through dimensions */
		for ( id = 0 ; id < ndim ; id++){
			/* clear X1 */
			path_scale(X1,X1,0,ndim+1);
			/* write dummy amounts */
			path_scale(X1, X, 1, ndim);
			X1 [ndim ] = (float)id*D1[ndim];
		    /* ealuate dS at this point, we are only allowing this to perturb  in first dimension*/
		    G [ ik*ndim + id ] = path_evaluate_array(X1, dPotential, N1, D1, O1 , ndim+1); 
		}
	}
	/* free arrays */
	free(X);
	free(X1);
	free( N1);
	free( O1);
	free (D1);
	
	return;
}

int path_find_tau_type( float V0, float V1, float V2)
	/*< finds the type of tau interpolation >*/
{
	int type;
	if ( V2 > V1 ){
		if ( V1 > V0){
			type = 1;
		}else{
			type = 0;
		}
	} else {
		if ( V1 < V0 ){
			type = 2;
		} else {
			type = 0;
		}
	}
	/* check which case for saddle point (0) */
	if ( type == 0 ){
		if ( V2 > V0 ){
			type = 3;
		}else{
			type = 4;
		}
	}
	return type;
}

float path_make_vimax ( float V0, float V1, float V2)
	/*< make vimax from the V values >*/
{
	return path_max( path_abs(V2-V1), path_abs(V1-V0));
}

float path_make_vimin ( float V0, float V1, float V2)
	/*< make vimin from the V values >*/
{
	return path_min( path_abs(V2-V1), path_abs(V1-V0));
}

void path_tau_type_to_tau_loc( float* TauL, float* Tau_pL, float* Tau_mL, float V0, float V1, float V2, int ndim)
	/*< creates tau local from tau plus and tau minus >*/
{
	float vimax = 0;
	float vimin = 0;
	/* find interpolation type */
	int type = path_find_tau_type(V0,V1,V2);
	if ( type == 1 ){
		/* Tau is taup */
		path_scale(TauL, Tau_pL, 1, ndim);
	}else{
		if ( type == 2 ){
		/* tau is taum*/
		    path_scale(TauL, Tau_mL, 1, ndim);
	    } else{
			/* determine Vimax and Vimin */
			vimax = path_make_vimax(V0, V1, V2);
			vimin = path_make_vimin(V0, V1, V2);
			if ( type == 3 ){
				path_combine(TauL, Tau_mL, vimin, Tau_pL, vimax, ndim);
			} else{
				if ( type == 4 ){
					path_combine(TauL, Tau_mL, vimax, Tau_pL, vimin, ndim);
				} else {
					sf_error("Invalid Tau Type during Runtime");
				}
			}
		}
	}
}

void path_create_tau_stable( float* Tau, float* Tau_p, float* Tau_m, float* V, int nknots, int ndim)
	/*< create tau in stable manner >*/
{
	/* looping index */
	int ik;
	/* local arrays*/
	float* Tau_pL = sf_floatalloc(ndim);
	float* Tau_mL = sf_floatalloc(ndim);
	float* TauL   = sf_floatalloc(ndim);
	/* zero out Tau */
	path_scale(Tau,Tau,0.,ndim*nknots);
	/* loop through knots */
	for ( ik = 0 ; ik < nknots ; ik++ ){
		/* get arrays */
		path_get_column(Tau_p, Tau_pL, ik, ndim);
		path_get_column(Tau_m, Tau_mL, ik, ndim);
		/* check end point cases */
		if ( ik == 0 ) {
			/* copy tau_p to tau */
			path_scale(TauL,Tau_pL,1,ndim);
		}else{
			if (ik == nknots-1){
				/* copy tau_m to tau*/
				path_scale(TauL,Tau_mL,1,ndim);
			} else {
				/* find tau local  */
				path_tau_type_to_tau_loc( TauL, Tau_pL, Tau_mL, V[ik-1], V[ik], V[ik+1], ndim);
			}
		}
		/* normalize local Tau */
		path_normalize(TauL,TauL,ndim);
		/* put Tau in array */
		path_put_column(Tau,TauL,ik,ndim);
	}
	free (Tau_pL);
	free (Tau_mL);
	free (TauL);
	return;
}

void path_create_tau(float* Tau, float* dR, int nknots, int ndim)
	/*<create the Tau approximation by averaging dR around point, this version is less stable>*/
{
	/* looping index */
	int i;
	/* local arrays */
	float* TauL = sf_floatalloc(ndim);
	float* dRo  = sf_floatalloc(ndim);
	float* dRf  = sf_floatalloc(ndim);
	/* zero out tau */
	path_scale( Tau, Tau, 0, nknots*ndim);
	for ( i = 1 ; i < nknots ; i++ ){
		/* get comuns */
		path_get_column(dR, dRo, i, ndim);
		path_get_column(dR, dRf, i+1, ndim);
		/* average */
		path_combine(TauL, dRo, 0.5, dRf, 0.5, ndim);
		/* put column back */
		path_put_column(Tau,TauL,i,ndim);
	}
	/* free arrays */
    	free (TauL);
		free ( dRo);
		free (dRf);
	return;
}

void path_create_tau_repulse(float* Tau, float* Spring, float* dR, int nknots, int ndim)
	/*<create the Tau approximation by averaging dR around point, also calculates precursor to spring force>*/
{
	/* looping index */
	int i;
	/* local arrays */
	float* TauL = sf_floatalloc(ndim);
	float* dRo  = sf_floatalloc(ndim);
	float* dRf  = sf_floatalloc(ndim);
	float* TauN = sf_floatalloc(ndim);
	float* SprL = sf_floatalloc(ndim);
	
	float dist ;
	
	/* zero out tau */
	path_scale( Tau, Tau, 0, nknots*ndim);
	for ( i = 0 ; i < nknots ; i++ ){
		/* get comuns */
		path_get_column(dR, dRo, i, ndim);
		path_get_column(dR, dRf, i+1, ndim);
		/* average */
		path_combine(TauL, dRo, 0.5, dRf, 0.5, ndim);
		/* put column back */
		path_put_column(Tau,TauL,i,ndim);
		/* normalize tau */
		path_normalize(TauN,Tau,ndim);
		/* get distance magnitude for spring */
		dist = path_norm(dRf,ndim) - path_norm( dRo, ndim);
		/* apply it to normalized path vector */
		path_scale(SprL,TauN,dist,ndim);
		/* put in array */
		path_put_column(Spring,SprL,i,ndim);
		
	}
	/* free arrays */
    	free (TauL);
		free ( dRo);
		free (dRf);
		free (TauN);
		free (SprL);
	return;
}


void path_orthagonalize_array(float* Orth, float* G, float* Tau, int knots, int ndim)
	/*< orthagonalize each component of G to Tau >*/
{
	/* local arrays */
	float* Gloc = sf_floatalloc(ndim);
	float* Orthloc = sf_floatalloc(ndim);
	float* Tauloc = sf_floatalloc(ndim);
	/* loop through knots */
	int ik;
	for (ik = 0 ; ik < knots ; ik++){
		/* get G and Tau */
		path_get_column(G,Gloc,ik,ndim);
		path_get_column(Tau,Tauloc,ik,ndim);
		/* orthagonalize */
		path_orthagonalize(Orthloc, Gloc, Tauloc, ndim);
		/* put it in Orth */
		path_put_column(Orth,Orthloc,ik,ndim);
	}

	free(Gloc);
	free(Orthloc);
	free(Tauloc);
	
	return;
}

void path_enforce_smoothness1(float* R, int order, int knots, int ndim)
	/*< enforces smoothness constraint on R , not currently functioning properly>*/
{
	/* loopoing index */
	int ik, io, this, lorder;
	float weight ;
	/* local arrays */
	float* PullX = sf_floatalloc(ndim);
	float* PuttX = sf_floatalloc(ndim);
	/* loop through */
	for ( ik = 0 ; ik < knots ; ik++){
		/* make sure we will be symmetric with our order */
		lorder = 0;
		for ( io = 0; io <= lorder ; io++){
			this = io + ik;
			if ( this > knots) continue;
			this = ik - io;
			if ( this < 0 ) continue;
			/* increment because we are OK! */
			lorder += 1;
		}
		/* zero out avg */
		path_scale(PuttX,PuttX,0,ndim);
		/* initialize counter */
		weight = 0 ;
		for ( io = -1*lorder; io <= lorder ; io++){
			if (io == 0 ) continue;
			this = io + ik;
			/* make sure this operation makes sense */
			if ( this < 0 || this > knots) continue;
			/* pull column */
			path_get_column(R,PullX,this,ndim);
			/* add to avg */
			path_combine(PuttX,PuttX,1,PullX,-1./path_abs((float)io),ndim);
			/* increment counter */
			weight += 1/path_abs((float)io);
		}
		
		/* scale by how many arrays we considered */
		path_get_column(R,PullX,ik,ndim);
		if (weight > 0 ){
			path_combine(PuttX,PuttX,1,PullX,-1*weight,ndim);
//			path_combine(PuttX,PuttX,1,PullX,1,ndim);
		} else{
			path_scale(PuttX, PullX, 1.0, ndim);
		}
		/* put it back */
		path_put_column( R, PuttX, ik, ndim);
	}
	free( PullX);
	free( PuttX);
	return;
}

void path_enforce_boundaries( float* R, int nknots, int* N, float* D, float* O, int ndim)
	/*< makes sure the R stays in bounds >*/
{
	/* looping index */
	int ik ;
	/* dimension index */
	int id ;
	/* do we need to put it back ? */
	bool writeswitch = false;
	/* local array */
	float* X = sf_floatalloc(ndim);
	/* create max array */
	float* M = sf_floatalloc(ndim);
	path_make_max( M, N, D, O, ndim);
	for ( ik = 0 ; ik < nknots ; ik++ ){
		/* get coordinate */
		path_get_column(R, X, ik, ndim);
		/* check to see if in bounds */
		for ( id = 0 ; id < ndim ; id++ ){
			/* are we too low? */
			if ( X[id] < O[ id]){
				writeswitch = true;
				X[id] = O[id];
			}
			/* are we too high ? */
			if ( X[id] > M[id]){
				writeswitch = true;
				X[id] = M[id];
			}
		}
		/* if we changed the array, put the column back */
		if (writeswitch){
			path_put_column( R, X, ik, ndim);
		}
	}
	free (X);
	free (M);
	return;
}

void path_enforce_boundaries_change( float* R, float* Change, float damp, int nknots, int* N, float* D, float* O, int ndim)
	/*< makes sure the R stays in bounds, reflects the change component and dampens >*/
{
	/* looping index */
	int ik ;
	/* dimension index */
	int id ;
	/* do we need to put it back ? */
	bool writeswitch = false;
	/* local array */
	float* X = sf_floatalloc(ndim);
	/* local change array */
	float* ChLoc = sf_floatalloc(ndim);
	/* create max array */
	float* M = sf_floatalloc(ndim);
	path_make_max( M, N, D, O, ndim);
	for ( ik = 0 ; ik < nknots ; ik++ ){
		/* get coordinate */
		path_get_column(R, X, ik, ndim);
		/* and get the change */
		path_get_column(Change,ChLoc,ik,ndim);
		/* check to see if in bounds */
		for ( id = 0 ; id < ndim ; id++ ){
			/* are we too low? */
			if ( X[id] < O[ id]){
				writeswitch = true;
				/* force to edge */
				X[id] = O[id];
				/* reflect and dampen change */
				ChLoc[id] = -1*damp*ChLoc[id];
			}
			/* are we too high ? */
			if ( X[id] > M[id]){
				writeswitch = true;
				/* force to edge */
				X[id] = M[id];
				/* reflect and dampen change */
				ChLoc[ id] = -1*damp*ChLoc[id];
			}
		}
		/* if we changed the array, put the column back */
		if (writeswitch){
			path_put_column( R, X, ik, ndim);
			path_put_column( Change, ChLoc, ik, ndim);
		}
	}
	free (X);
	free (M);
	free (ChLoc);
	return;
}

void path_enforce_function( float* R, float eps, int knots, int ndim)
	/*< enforce requirement that R(t) >*/
{
	/* looping index */
	int ik;
	/* local indexes */
	float* X0 = sf_floatalloc(ndim);
	float* X1 = sf_floatalloc(ndim);
	for ( ik = 1 ; ik < knots ; ik++ ){
		/* get local arrays */
		path_get_column(R, X0, ik-1, ndim);
		path_get_column(R, X1, ik, ndim);
		if ( X0[0] > X1[0]){
			X1[0] = X0[0] + eps;
		}
		/* put it back */
		path_put_column(R,X1,ik,ndim);
	}
	free(X0);
	free(X1);
	return;
}

void path_unkink_force( float* Unkink, float* Tau_p, float* Tau_m, int nknots, int ndim)
	/*< a force that resists bends in the chain >*/
{
	/* looping index */
	int ik;
	/* local array */
	float* Tau_pL = sf_floatalloc(ndim);
	float* Tau_mL = sf_floatalloc(ndim);
	float* Tau_dif = sf_floatalloc(ndim);
	float* Unkink_L = sf_floatalloc(ndim);
	/* zero out unkink */
	path_scale (Unkink, Unkink, 0.,nknots*ndim);
	/* loop through knots */
	for ( ik = 0 ; ik < nknots ; ik++){
		/* get columns */
		path_get_column( Tau_p, Tau_pL, ik, ndim);
		path_get_column( Tau_m, Tau_mL, ik, ndim);
		path_combine( Tau_dif, Tau_pL, 1, Tau_mL, -1, ndim);
		path_put_column( Unkink, Tau_dif, ik, ndim);
	}
	/* free arrays */
	free( Tau_pL);
	free( Tau_mL);
	free( Unkink_L);
	
	return;
}


void path_attractive_force(float* Attract, float* G, float* Tau_p, float* Tau_m, int nknots, int ndim)
	/*< create the attractive force felt by the path >*/
{
	/* looping index */
	int ik;
	/* local arrays */
	float* Gloc   = sf_floatalloc(ndim);
	float* Tau_pL = sf_floatalloc(ndim);
	float* Tau_mL = sf_floatalloc(ndim);
	float* At_loc   = sf_floatalloc(ndim);
	/* zero out attract */
	path_scale( Attract, Attract, 0, nknots*ndim);
	/* distance scaling of gradient */
	float dist ; 
	/* loop through knots */
	for ( ik = 0 ; ik < nknots ; ik++ ){
		/* get local arrays */
		path_get_column(G, Gloc, ik, ndim);
		path_get_column(Tau_p, Tau_pL, ik, ndim);
		path_get_column(Tau_m, Tau_mL, ik, ndim);
		/* calculate distance of string to adjacent nodes */
		dist = path_norm ( Tau_pL, ndim) + path_norm ( Tau_mL, ndim);
		/* scale gradient by distance */
		path_scale( At_loc, Gloc, dist, ndim);
		/* force the first and last entries to have zero derivative in x1 direction */
		if ( ik == 0 ){
			At_loc[ 0] = 0;
		}
		if ( ik == nknots-1){
			At_loc[ 0] = 0;
		}
		/* put the force in the global array */
		path_put_column( Attract, At_loc, ik, ndim);
	}
	
	free( Gloc);
	free( Tau_pL);
	free(Tau_mL);
	free(At_loc);
	return;
}

void path_repulsive_force(float* Repulse, float* V, float* Tau_p, float* Tau_m, int nknots, int ndim)
	/*< get the repulsive force between nodes in the path >*/
{
	/* looping index */
	int ik;
	/* local arrays */
	float* Tau_pL = sf_floatalloc(ndim);
	float* Tau_mL = sf_floatalloc(ndim);
	float* Rep_L  = sf_floatalloc(ndim);
	/* scalings for directions */
	float scale_p, scale_m;
	/* zero out repulsion */
	path_scale(Repulse, Repulse, 0, nknots*ndim);
	/* loop through knots */
	for ( ik = 0 ; ik < nknots ; ik++ ){
		/* get local arrays */
		path_get_column(Tau_p, Tau_pL, ik, ndim);
		path_get_column(Tau_m, Tau_mL, ik, ndim);
		/* normalize */
		path_normalize ( Tau_p, Tau_p, (long)ndim);
		path_normalize ( Tau_m, Tau_p, (long)ndim);
		/* get scaling factors */
		if ( ik > 0 ){
			scale_m = V[ ik] - V[ ik-1];
		} else {
			scale_m = 0.;
		}
		if ( ik < nknots - 1 ){
			scale_p = V[ ik+1] - V[ik];
		} else {
			scale_p = 0.;
		}
		/* combine */
		path_combine( Rep_L, Tau_m, -1*scale_m, Tau_p, scale_p, ndim);
		/* put in global array */
		path_put_column(Repulse, Rep_L, ik, ndim);
	}
	free (Rep_L);
	free (Tau_pL);
	free (Tau_mL);
	return;
}

void path_spring_force( float* Spring, float* Tau, float* Tau_p, float* Tau_m, int nknots, int ndim)
	/*< create spring forces >*/
{
	/* looping index */
	int ik;
	/* local arrays */
	float* Tau_pL = sf_floatalloc(ndim);
	float* Tau_mL = sf_floatalloc(ndim);
	float* Tau_L  = sf_floatalloc(ndim);
	float* Spr_L  = sf_floatalloc(ndim);
	/* distance scales */
	float dist_m, dist_p, dist;
	
	/* loop through knots */
	for ( ik = 0 ; ik < nknots; ik++ ){
		/* get local arrays */
		path_get_column(Tau_p, Tau_pL, ik, ndim);
		path_get_column(Tau_m, Tau_mL, ik, ndim);
		path_get_column(  Tau, Tau_L , ik, ndim);
		if ( ik > 0 ){
			dist_m = path_norm(Tau_mL, ndim);
		} else {
			dist_m = 0.;
		}
		if ( ik < nknots - 1 ){
			dist_p = path_norm(Tau_pL, ndim);
		} else {
			dist_p = 0.;
		}
		/* subtract */
		dist = (dist_p - dist_m);
		/* create spring force */
		path_scale(Spr_L, Tau_L, dist, ndim);
		/* put array */
		path_put_column(Spring, Spr_L, ik, ndim);
	}
	free(Tau_pL);
	free(Tau_mL);
	free(Tau_L);
	free(Spr_L);
	return;
}


void path_reinterpolate1(float* Interpolated, float* R, int* N, float* D, float* O, int nknots, int ndim)
	/*< interpolates our path knots to the sampling rate >*/
{
	/* looping index */
	int i1;
	/* local arrays, first and second knot */
	float* R0 = sf_floatalloc(ndim);
	float* R1 = sf_floatalloc(ndim);
	/* zero out interpolated array */
	path_scale(Interpolated,Interpolated,0.,N[0]);
	
	float x1, rx1dst;
	/* which knot are we on?*/
	int knot = 0;
	/* get initial knots */
	path_get_column(R,R0,knot,ndim);
	path_get_column(R,R1,knot+1,ndim);
	rx1dst = R1[0] - R0[0];
	/* set first point */
	/* loop through samples */
	for (i1 = 0 ; i1 < N[0] ; i1++ ){
		/* where are we ?*/
		x1 = ((float)i1) * D[0] + O[0];
		/* do we need new knots ?*/
		if ( x1 > R1[0] ){
			/* get next knot */
			knot += 1;
			path_get_column(R,R0,knot,ndim);
			if (knot+1 < nknots){ 
				path_get_column(R,R1,knot+1,ndim);
			}else{ break; }
			if ( x1 > R1[0]){
				i1 -= 1;
				continue;
			}
			/* distance between x components */
			rx1dst = R1[0] - R0[0];
		}
		/* interpolate */
		Interpolated [ i1] = ((x1-R0[0]) * R1[1] + (R1[0]-x1)* R0[1])/rx1dst;
	}
	/* free arrays */
	free ( R0);
	free ( R1);
	return;
}
	

