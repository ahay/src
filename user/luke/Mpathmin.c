/* One dimensional program for iteratively determining minimal traveltime path for
   through a cost function input file has as the first dimension a parameter that can vary in the path, 
   second dimension "time", or a coordinate that the output path, R(t), will be a function of. 
   Higher dimensions represent additional panels to calculate the path through.
   To use this program for picking a path maximizing an objective function, first transform the 
   objective function \alpha by exp(-1*\alpha) and feed the resulting output to this program.
   Path integral minimization is in the manner of "Path Optimization with Application to Tunneling",
   Dorothea M. Einarsdottir, Andri Arnaldsson, Finnbogi Oskarsson, and Hannes Jonsson, 2012.
*/

#include <rsf.h>
#include <stdlib.h>
#include <time.h>
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
	
    /* Get file sampling info */
	int n1, n2, nmid;
	int id, ik;
	float d1, d2;
	float o1, o2;
    if (SF_FLOAT != sf_gettype(_in)) sf_error("Need float input");
    fildim = sf_filedims (_in,NFile);
    if (fildim < 2) sf_error("Need at least two dimensions");
	
	/* copy dimension info from file */
	n1 = NFile[0]; 
	n2 = NFile[1]; 
	
	/* put all of the midpoints (nx, ny, ...) on n3 if we are calculating path on more than one panel*/
	nmid = 1;
	for (id=2; id < fildim; id++) nmid *= NFile[id];
	/* sampling info */
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
	
	/* shift dimensions in output (2 dim panel -> 1 dim path)*/
	sf_unshiftdim(_in,_out,2);	
	/* label info */
    if (NULL != (label = sf_histstring(_in,"label2")))
		sf_putstring(_out,"label",label);
    if (NULL != (label = sf_histstring(_in,"unit2")))
		sf_putstring(_out,"unit",label);
	
	/* get path optimization parameters */
	
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
	
	int niter ; 
	if (!sf_getint("niter",&niter))   niter = 10;
	/* number of iterations */
	
	float damp;
	if (!sf_getfloat("damp",&damp))   damp = .5;
	/* if the path goes out of bounds, we reflect and dampen the rate of change by this much */	
	
	float shove;
	if (!sf_getfloat("shove",&shove))   shove = 1000;
	/* size of initial random lateral shove */	
			
	/* allocate sampling arrays, only using ndim=2 for 2d picking panels here*/
	/* number of samples per dimension*/
	int* N = sf_intalloc(ndim);
	N[0] = n1;
	N[1] = n2;
	
	/* sampling increment */
	float* D = sf_floatalloc(ndim);
	D[0] = d1;
	D[1] = d2;
	
	/* axis origin */
	float* O = sf_floatalloc(ndim);
	O[0] = o1;
	O[1] = o2;
	
	float aniso1;
	if (!sf_getfloat("aniso1",&aniso1)) aniso1=D[1]/D[0];
	/* anisotropy of 2nd axis relative to first   */
	
	/* input panel gradient calculation parameters */
	int dorder ;
	if (!sf_getint("dorder",&dorder)) dorder=6;
	/* derivative order (stencil size) for gradient calculation*/

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
	
	/* scaling gradient by 1/D[i]? Otherwise assumes sampling interval D[i] = 1*/
	bool scale = true;
	
	/* orth force magnitude */
	float orthnorm;
	/* spring force magnitude */
	float sprnorm;
	/* kink force magnitude */
	float kinknorm;
	/* size of a selection panel */
	long panelsize = path_size(N,ndim);
	/* allocate S, the input panel we pick on */
	float* Sfunc = sf_floatalloc(panelsize);
	/* and array for its gradient */
	float* dSfunc = sf_floatalloc(panelsize*ndim);
	/* interpolated path to be output */
	float* Interpolated = sf_floatalloc(N[0]);
	/* allocate initial path, R */
	float* R = sf_floatalloc(knots*ndim);
	/* allocate potential at each knot in R*/
	float* V = sf_floatalloc(knots);
	/* allocate gradient of S at each knot in R*/
	float* G = sf_floatalloc(knots*ndim);
	/* allocate attractive force at each knot in R*/
	float* Attract = sf_floatalloc(knots*ndim);
	/* allocate repulsive force at each knot in R */
	float* Repulse = sf_floatalloc(knots*ndim);
	/* allocate combined forces at each knot in R */
	float* Force   = sf_floatalloc(knots*ndim);
	/* allocate force orthagonal to path R at each knot in R */
	float* Orth    = sf_floatalloc(knots*ndim);
	/* allocate unkinking force at each knot in R*/
	float* Unkink  = sf_floatalloc(knots*ndim);
	/* allocate spring force pushing knots apart at each knot in R */
	float* Spring  = sf_floatalloc(knots*ndim);
	/* allocate spring force scaled by stiffness k at each knot in R */
	float* SpringK  = sf_floatalloc(knots*ndim);
	/* allocate array for change of gradient at each knot in R */
	float* Change = sf_floatalloc(knots*ndim);
	/* allocate array the update to R at each knot */
	float* Update = sf_floatalloc(knots*ndim);
	/* g scaled update */
	float* UpdateG = sf_floatalloc(knots*ndim);
	
	/* allocate path parallel vector and its precursors for determining orth and component in manner of: */
	/* Improved tangent estimate in the nudged elastic band method for finding minimum energy 
	   paths and saddle points, Graeme Henkelmanan and Hannes Jonsson, 2000 */
	float* Tau = sf_floatalloc(knots*ndim);
	/* tau plus */
	float* Tau_p = sf_floatalloc(knots*ndim);
	/* tau minus */
	float* Tau_m = sf_floatalloc(knots*ndim);
	/* seed random number generator */
	srand(time(0)*time(0));	
	float randum ;
	/* loop through midpoints */
	int imid, iter;
	for ( imid = 0 ; imid < nmid ; imid ++ ){
		sf_warning("Learning path %i of %i;",imid+1,nmid);
		/* get intial guess for path R */
		path_first_guess1(R, knots, N, D, O, ndim);
		/* make first update vector*/
		randum = ((float)rand()/(.5*(float)RAND_MAX)-1)*shove;
		path_first_update1(Update,randum,knots,ndim);
		/* make sure we are in bounds */
		path_enforce_boundaries( R, knots, N, D, O, ndim);	
		/* zero out S */
		path_scale( Sfunc, Sfunc, 0., panelsize);
		/* zero out the grad of S */
		path_scale( dSfunc, dSfunc, 0., panelsize*ndim);	
		/* zero out V */
		path_scale( V, V, 0., knots);
		/* zero out G */
		path_scale(G, G, 0., knots*ndim);
		/* zero out attractive force */
		path_scale( Attract, Attract, 0., knots*ndim);
		/* zero out Repulse force */
		path_scale( Repulse, Repulse, 0, knots*ndim);
		/* zero out total force */
		path_scale(Force, Force, 0., knots*ndim);
		/* zero out force orthagonal to path */
		path_scale(Orth, Orth, 0., knots*ndim);
		/* zero out unkinking force */
		path_scale(Unkink, Unkink, 0., knots*ndim);
		/* zero outo spring force */
		path_scale(Spring, Spring, 0., knots*ndim);
		/* zero out spring times K */
		path_scale(SpringK, SpringK, 0., knots*ndim);
		/* zero out change array */
		path_scale(Change, Change, 0., knots*ndim);
		/* zero out G update  */
		path_scale(UpdateG, UpdateG, 0., knots*ndim);
		/* zero out interpolated path for output*/
		path_scale(Interpolated,Interpolated, 0.,N[0]);
		/* read S, the panel we are traveling through, from file */
		sf_floatread(Sfunc,panelsize,_in);
		/* compute gradient */
		path_gradient( Sfunc, dSfunc,  N, D, O, dorder, slen, nsmooth, ndim, scale);
		/* evaluate potential (S) at knot points */
		path_evaluate_potental(V, R, knots, Sfunc, N, D, O, ndim);
		/* evaluate gradient of S at knot points */
		path_evaluate_gradient(G, R, knots, dSfunc, N, D, O, ndim);
		
		/* iterate to learn path */
		for ( iter = 0 ; iter < niter ; iter++ ){
			/* get tau plus and tau minus (un-normalized), 
			   precursor to actual tau, the tangent */
			path_create_tau_plus_minus(Tau_p, Tau_m, R, knots, ndim); 
			/* combine data into tangent curve Tau */
			path_create_tau_stable( Tau, Tau_p, Tau_m, V, knots, ndim);
			/* apply anisotropy to make a preferred direction.  
			   Higher anisotroy makes path more willing to bend "into" picking coordinate
			   and thus become more serpentine */
			for ( ik = 0 ; ik < knots ; ik++ ){
				Tau_p[ik*ndim +1 ] *= aniso1;
				Tau_m[ik*ndim +1 ] *= aniso1;
			}
			/* unkinking force that pushes back against bends in path */
			path_unkink_force( Unkink, Tau_p, Tau_m, knots, ndim);
			/* get the attractive force pulling knots together trying to minimize path cost 
			   by "pulling" the path into lows of S along gradient */
			path_attractive_force( Attract, G, Tau_p, Tau_m, knots, ndim);
			/* get the repulsive force trying to keep the integral over path length smaller */
			path_repulsive_force( Repulse, V, Tau_p, Tau_m, knots, ndim);
			/* combine forces */
			path_combine( Force, Attract, -1, Repulse, 1, knots*ndim);
			/* orthagonalize force  so this component only works perpendicular to path*/
			path_orthagonalize_array( Orth, Force, Tau, knots, ndim);
			/* force first and last knot's force to be 0 */
			for ( id = 0 ; id < ndim; id++){
				Orth[ id] = 0.;
				Orth[ (knots-1)*ndim +id ] = 0.;
			}
			/* calculate magnitude of orthagonal force */
			orthnorm = path_norm(Orth, ndim*knots);
			/* create spring force pushing knots apart */
			path_spring_force( Spring, Tau, Tau_p, Tau_m, knots, ndim);
			/* compute magnitudes of spring and kink forces */
			sprnorm = path_norm(Spring, ndim*knots);
			kinknorm = path_norm(Unkink, ndim*knots);
			/* is the spring force magnitude nonzero (does division make sense)? */
			if (sprnorm > 0){
				k2 = k*orthnorm/sprnorm;
			} else {
				k2 = 0;
			}
			/* make the unkinking force mag proportional to the spring force magnitude */
			if ( kinknorm > 0){
				kink2 = kink * k2;
			}
			/* combine, scale by k */
			path_combine(SpringK, Spring, k2, Unkink, kink2, knots*ndim);
			/* combine forces into change of gradient (for the Path! not S!!) */
			path_combine( Change, SpringK, 1, Orth, 1, knots*ndim);
			/* enforce edges, cant pull the initial or last knots off of t=0, t=tf */
			Change[0] = 0;
			Change[(knots-1)*ndim] = 0;
			/* update the gradient of how R evolves with this iteration's calculated change */
			/* scaled by the learning rate */
			path_combine( Update, Update, 1-lr, Change, lr, knots*ndim);
			/* scale update by g before applying it to the knots */
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
			/* evaluate potential at knot points */
			path_evaluate_potental(V, R, knots, Sfunc, N, D, O, ndim);
			/* evaluate gradient at knot points */
			path_evaluate_gradient(G, R, knots, dSfunc, N, D, O, ndim);
		}
		/* enforce that the path must be function of t before we interpolate */
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

