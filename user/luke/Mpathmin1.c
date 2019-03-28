/* One dimensional path minimization for optimization input file has first coordinate parameter, second coordinate time */

#include <rsf.h>
#include "path.h"

int main (int argc, char* argv[])
{
	/* how many dimensions ?*/
	int ndim = 2;
	
	/* declare files */
    sf_file _in,  _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* initialize files */
	_in  = sf_input ( "in");
	_out = sf_output("out");
	
	int n1, n2;
	float d1, d2;
	float o1, o2;
	
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
	
	float k;
	float k2;
	if (!sf_getfloat("k",&k))   k = 1;
	/* stiffness relative to attraction*/
	
	float kink;
	float kink2;
	if (!sf_getfloat("kink",&kink)) kink=1;
	/* resistance to kinks  */
	
	float lr;
	if (!sf_getfloat("lr",&lr))   lr = .3;
	/* learning rate */
	if ( lr < 0 ){ 
		lr = .5;
		sf_warning("No negative Learning Rates, changed to 0.5");
	}
	if ( lr > 1 ){
		sf_warning(" Too high a learning rate, changed to 1");
		lr = 1;
	}

	float g;
	if (!sf_getfloat("g",&g))   g = .1;
	/* scaling the gradient by how much */	
	
	int knots;
    if (!sf_getint("knots",&knots))   knots = 11;
    /* number of knots */
	
	/* number of iterations */
	int niter ; 
	if (!sf_getint("niter",&niter))   niter = 10;
	
	/* initialize output */
	/* for writing out knots */
/*	sf_putint   (_out,"n1",ndim); 
	sf_putfloat (_out,"d1",1);
	sf_putfloat (_out,"o1",0);	
	sf_putint   (_out,"n2",knots); 
	sf_putfloat (_out,"d2",1.);
	sf_putfloat (_out,"o2",0.);
	*/
	/* for writing out interpolated result */
	sf_putint   (_out,"n1",n1); 
    sf_putfloat (_out,"d1",d1);
    sf_putfloat (_out,"o1",o1);	
	sf_putint   (_out,"n2",1 ); 
	sf_putfloat (_out,"d2",1.);
	sf_putfloat (_out,"o2",0.);
	
	/* for writing out gradient */
/*	sf_putint   (_out,"n1",n1); 
    sf_putfloat (_out,"d1",d1);
    sf_putfloat (_out,"o1",o1);	
	sf_putint   (_out,"n2",n2); 
    sf_putfloat (_out,"d2",d2);
    sf_putfloat (_out,"o2",o2);	
	sf_putint   (_out,"n3",ndim); 
    sf_putfloat (_out,"d3",1);
    sf_putfloat (_out,"o3",0);	
*/				
	/* allocate sampling arrays */
	int* N = sf_intalloc(ndim);
	N[0] = n1;
	N[1] = n2;
	
	float* D = sf_floatalloc(ndim);
	D[0] = d1;
	D[1] = d2;
	
	float* O = sf_floatalloc(ndim);
	O[0] = o1;
	O[1] = o2;
	
	
	float aniso1;
	if (!sf_getfloat("aniso1",&aniso1)) aniso1=D[1]/D[0];
	/* aniso of 2nd axis relative to first   */
	
	/* general purpose dimension index */
	int id;
	/* and knots */
	int ik;
	
	int dorder ;
	if (!sf_getint("dorder",&dorder)) dorder=6;
	/* derivative order */


	int slen;
	if (!sf_getint("srad",&slen)) slen=2;
	/* smoothing radius for gradient */
	
	int nsmooth;
	if (!sf_getint("nsmooth",&nsmooth)) nsmooth=1;
	/* number of gradient smoothings  */
	
	/* orth size */
	float orthnorm;
	/* spring size */
	float sprnorm;
	/* kink size */
	float kinknorm;
	/* size of a selection panel */
	long panelsize = path_size(N,ndim);
	/* allocate S */
	float* Sfunc = sf_floatalloc(panelsize);
	/* and array for its derivative, only 1d for now */
	float* dSfunc = sf_floatalloc(panelsize*ndim);
	/* read S from file */
	sf_floatread(Sfunc,panelsize,_in);
	/* compute gradient */
	/* scaling by 1/D[i]? */
	bool scale = true;
	path_gradient( Sfunc, dSfunc,  N, D, O, dorder, slen, nsmooth, ndim, scale);
	/* generate initial path */
	float* R = sf_floatalloc(knots*ndim);
	/* take first guess */
	path_first_guess1(R, knots, N, D, O, ndim);
	/* make sure we are in bounds */
	path_enforce_boundaries( R, knots, N, D, O, ndim);
	/* potential */
	float* V = sf_floatalloc(knots);
	/* grdient */
	float* G = sf_floatalloc(knots*ndim);
	/* attractive force */
	float* Attract = sf_floatalloc(knots*ndim);
	/* repulsive force */
	float* Repulse = sf_floatalloc(knots*ndim);
	/* combined force */
	float* Force   = sf_floatalloc(knots*ndim);
	/* orth force */
	float* Orth    = sf_floatalloc(knots*ndim);
	/* unkinking force */
	float* Unkink  = sf_floatalloc(knots*ndim);
	/* spring force */
	float* Spring  = sf_floatalloc(knots*ndim);
	/* spring force scaled by k*/
	float* SpringK  = sf_floatalloc(knots*ndim);
	/* change of gradient */
	float* Change = sf_floatalloc(knots*ndim);
	/* zero out just to be sure*/
	path_scale(Change, Change, 0., knots*ndim);
	/* the update for R */
	float* Update = sf_floatalloc(knots*ndim);
	/* zero out just to be sure */
	path_scale(Update, Update, 0., knots*ndim);
	/* g scaled update */
	float* UpdateG = sf_floatalloc(knots*ndim);
	/* zero out just to be sure */
	path_scale(UpdateG, UpdateG, 0., knots*ndim);
	/* line parallel vector and its precursors */
	float* Tau = sf_floatalloc(knots*ndim);
	/* tau plus */
	float* Tau_p = sf_floatalloc(knots*ndim);
	/* tau minus */
	float* Tau_m = sf_floatalloc(knots*ndim);	

	int iter ;

	for ( iter = 0 ; iter < niter ; iter++ ){
		/* evaluate potential at knot points */
		path_evaluate_potental(V, R, knots, Sfunc, N, D, O, ndim);
		/* evaluate gradient at knot points */
		path_evaluate_gradient(G, R, knots, dSfunc, N, D, O, ndim);
		/* get tau plus and tau minus (un-normalized), precursor to actual tau */
		path_create_tau_plus_minus(Tau_p, Tau_m, R, knots, ndim); // this is where a problem is!!!!!!!!!!!!
		/* combine data into tangent curve Tau */
		path_create_tau_stable( Tau, Tau_p, Tau_m, V, knots, ndim);
		/* apply anisotropy */
		for ( ik = 0 ; ik < knots ; ik++ ){
			Tau_p[ik*ndim +1 ] *= aniso1;
			Tau_m[ik*ndim +1 ] *= aniso1;
		}
		/* get the attractive force */
		path_attractive_force( Attract, G, Tau_p, Tau_m, knots, ndim);
		/* get the repulsive force */
		path_repulsive_force( Repulse, V, Tau_p, Tau_m, knots, ndim);
		/* combine forces */
		path_combine( Force, Attract, -1, Repulse, 1, knots*ndim);
		/* orthagonalize force */
		path_orthagonalize_array( Orth, Force, Tau, knots, ndim);
		/* force first and last to be 0 */
		for ( id = 0 ; id < ndim; id++){
			Orth[ id] = 0.;
			Orth[ (knots-1)*ndim +id ] = 0.;
		}
		/* calculate size of orthagonal force */
		orthnorm = path_norm(Orth, ndim*knots);
		/* create spring force */
		path_spring_force( Spring, Tau, Tau_p, Tau_m, knots, ndim);
		/* unkinking force */
		path_unkink_force( Unkink, Tau_p, Tau_m, knots, ndim);
		sprnorm = path_norm(Spring, ndim*knots);
		kinknorm = path_norm(Unkink, ndim*knots);
		if (sprnorm > 0){
			k2 = k*orthnorm/sprnorm;
		} else {
			k2 = 0;
		}
		if ( kinknorm > 0){
			kink2 = kink * k2;
		}
		/* combine, scale by k */
		path_combine(SpringK, Spring, k2, Unkink, kink2, knots*ndim);
		/* combine forces into change of gradient */
		path_combine( Change, SpringK, 1, Orth, 1, knots*ndim);
		/* enforce edges */
		Change[0] = 0;
		Change[(knots-1)*ndim] = 0;
		/* update the gradient with the change */
		path_combine( Update, Update, 1-lr, Change, lr, knots*ndim);
		/* scale by g */
		path_scale(UpdateG, Update, g, knots*ndim);
		/* and now lets change our knots */
		path_combine( R, R, 1., UpdateG, 1, knots*ndim);

		/* keep the end points on the edge */
		R [ 0] = O[ 0];
		R [ (knots-1)*(ndim)] = O[ 0] + D[ 0]*( (float) N[ 0] -1 );
		/* make sure we are in bounds */
		path_enforce_boundaries( R, knots, N, D, O, ndim);
	}
	/* enforce function before we interpolate */
	path_enforce_function( R, D[0], knots, ndim);
 //   path_enforce_boundaries( R, knots, N, D, O, ndim);
	float* Interpolated = sf_floatalloc(N[0]);
	/* and interpolate */
	path_reinterpolate1(Interpolated, R, N, D, O, knots, ndim);
	/* write interpolated to disk */
	sf_floatwrite(Interpolated,N[0],_out);
	
	/* write out result */
//	sf_floatwrite(R,knots*ndim,_out);
	//		sf_floatwrite(dSfunc,panelsize*ndim,_out);
	/* close file */
	sf_fileclose(_out);
	/* free arrays */
	free (Interpolated);
	free (Tau);
	free (Tau_p);
	free (Tau_m);
	free (G);
	free (V);
	free ( Sfunc);
	free (dSfunc);
	free (R);
	free( N);
	free( D);
	free( O);
	free (Attract);
	free (Repulse);
	free (Orth);
	free (Spring);
	free (SpringK);
	free (UpdateG);
	/* done here */
	exit(0);
}

