/* One dimensional path minimization for optimization input file has first coordinate parameter, second coordinate time */

#include <rsf.h>
#include "path.h"

int main (int argc, char* argv[])
{
	/* how many dimensions are we picking over ?*/
	int ndim = 2;
	
	/* dimensions of file size */
	int fildim, NFile[SF_MAX_DIM];
	char *label;
	
	/* declare files */
    sf_file _in,  _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* initialize files */
	_in  = sf_input ( "in");
	_out = sf_output("out");
	
	int n1, n2, nmid;
	int id, ik;
	float d1, d2;
	float o1, o2;
	
    /* Get sampling info */

	
    if (SF_FLOAT != sf_gettype(_in)) sf_error("Need float input");
    fildim = sf_filedims (_in,NFile);
    if (fildim < 2) sf_error("Need at least two dimensions");
	
	/* copy to variables */
	n1 = NFile[0]; 
	n2 = NFile[1]; 
	/* put all of the midpoints on n3 */
	nmid = 1;
	for (id=2; id < fildim; id++) nmid *= NFile[id];

	/* sampling info */
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
	
	/* shift dimensions in output*/
	sf_unshiftdim(_in,_out,2);
		
	/* label info */
    if (NULL != (label = sf_histstring(_in,"label2")))
		sf_putstring(_out,"label",label);
    if (NULL != (label = sf_histstring(_in,"unit2")))
		sf_putstring(_out,"unit",label);
	
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
	
	float damp;
	if (!sf_getfloat("damp",&damp))   damp = .5;
	/* if the path goes out of bounds, we reflect and dampen the rate of change by this much */	
			
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
	
	int dorder ;
	if (!sf_getint("dorder",&dorder)) dorder=6;
	/* derivative order */


	int slen;
	if (!sf_getint("srad",&slen)) slen=2;
	/* smoothing radius for gradient */
	
	int nsmooth;
	if (!sf_getint("nsmooth",&nsmooth)) nsmooth=1;
	/* number of gradient smoothings  */
	
	/* termination information */
	float update_size, change_size, termU, termC;
	
	float eps;
	if (!sf_getfloat("eps",&eps)) eps=0.;
	/* if the change and gradient are simultaneously lower than this, terminate  early */
	/* set both termination parameters to epsilon */
	termU = eps;
	termC = eps;
	
	/* scaling gradient by 1/D[i]? */
	bool scale = true;
	
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
	/* interpolated path */
	float* Interpolated = sf_floatalloc(N[0]);
	/* generate initial path */
	float* R = sf_floatalloc(knots*ndim);
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
	/* the update for R */
	float* Update = sf_floatalloc(knots*ndim);
	/* g scaled update */
	float* UpdateG = sf_floatalloc(knots*ndim);
	/* line parallel vector and its precursors */
	float* Tau = sf_floatalloc(knots*ndim);
	/* tau plus */
	float* Tau_p = sf_floatalloc(knots*ndim);
	/* tau minus */
	float* Tau_m = sf_floatalloc(knots*ndim);	
	/* loop through midpoins */
	int imid, iter;
	for ( imid = 0 ; imid < nmid ; imid ++ ){
		sf_warning("Learning path %i of %i;",imid+1,nmid);
		/* get intial guess for R */
		path_first_guess1(R, knots, N, D, O, ndim);
		/* make sure we are in bounds */
		path_enforce_boundaries( R, knots, N, D, O, ndim);	
		/* zero out S */
		path_scale( Sfunc, Sfunc, 0., panelsize);
		/* zero out the grad S */
		path_scale( dSfunc, dSfunc, 0., panelsize*ndim);	
		/* zero out V */
		path_scale( V, V, 0., knots);
		/* zero out G */
		path_scale(G, G, 0., knots*ndim);
		/* zero out attract */
		path_scale( Attract, Attract, 0., knots*ndim);
		/* zero out Repulse */
		path_scale( Repulse, Repulse, 0, knots*ndim);
		/* zero out force */
		path_scale(Force, Force, 0., knots*ndim);
		/* zero out orth */
		path_scale(Orth, Orth, 0., knots*ndim);
		/* zero out unkink */
		path_scale(Unkink, Unkink, 0., knots*ndim);
		/* zero outo spring */
		path_scale(Spring, Spring, 0., knots*ndim);
		/* zero out spring K */
		path_scale(SpringK, SpringK, 0., knots*ndim);
		/* zero out Change */
		path_scale(Change, Change, 0., knots*ndim);
		/* zero out Update  */
		path_scale(Update, Update, 0., knots*ndim);
		/* zero out G  */
		path_scale(UpdateG, UpdateG, 0., knots*ndim);
		/* zero out interpolated path */
		path_scale(Interpolated,Interpolated, 0.,N[0]);
		/* read S from file */
		sf_floatread(Sfunc,panelsize,_in);
		/* compute gradient */
		path_gradient( Sfunc, dSfunc,  N, D, O, dorder, slen, nsmooth, ndim, scale);
		/* iterate to learn path */
		for ( iter = 0 ; iter < niter ; iter++ ){
			/* evaluate potential at knot points */
			path_evaluate_potental(V, R, knots, Sfunc, N, D, O, ndim);
			/* evaluate gradient at knot points */
			path_evaluate_gradient(G, R, knots, dSfunc, N, D, O, ndim);
			/* get tau plus and tau minus (un-normalized), precursor to actual tau */
			path_create_tau_plus_minus(Tau_p, Tau_m, R, knots, ndim); 
			/* combine data into tangent curve Tau */
			path_create_tau_stable( Tau, Tau_p, Tau_m, V, knots, ndim);
			/* apply anisotropy */
			for ( ik = 0 ; ik < knots ; ik++ ){
				Tau_p[ik*ndim +1 ] *= aniso1;
				Tau_m[ik*ndim +1 ] *= aniso1;
			}
			/* unkinking force */
			path_unkink_force( Unkink, Tau_p, Tau_m, knots, ndim);
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
			/* make sure we are in bounds, reflect and dampen the update if at edge */
			path_enforce_boundaries_change( R, Update, damp, knots, N, D, O, ndim);
		
			/* determine if we are terminating because of convergence */
			update_size = path_norm(Update,knots*ndim);
			change_size = path_norm(Change,knots*ndim);
			if ( update_size < termU && change_size < termC ){
				/* terminate */
				sf_warning("Path %i learned in %i iterations",imid+1,iter+1);
				break;
			}
		}
		/* enforce function before we interpolate */
		path_enforce_function( R, D[0], knots, ndim);
		/* and interpolate */
		path_reinterpolate1(Interpolated, R, N, D, O, knots, ndim);
		/* write interpolated to disk */
		sf_floatwrite(Interpolated,N[0],_out);
	}
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
	free (label);
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

