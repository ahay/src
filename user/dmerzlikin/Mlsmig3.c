/* Least-Squares 3D Path-Summation Integral, Azimuthal Plane-Wave Destruction and Kirchhoff Modeling/Migration Chain of Operators*/
#include <rsf.h>
#include <math.h>
#include "testmig3.h"
#include "anisodiffuse.h"
#include "thr.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    int nt, nt2, nx,ny, i1, i2, i3, n12,nt12, i, j, initer, oniter, niter, repeat, dsnaps;
    bool adj, sm, domod, dopi, doanisodiff, dothr, thrflag, doomp, snaps, ch=false;
    float dt, dt2, dx,dy, ot, ot2, ox,oy, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr, thr;
    float *data, *modl, *modl0;
    sf_file inp, out;
    /* PWD parameters */
    int nw, nj1, nj2, apt;
    float *pp1, *pp2, *az;
    float anisoeps, *vx, *vy;
    sf_file dip, azin;
    /* Kirchhoff params */
    bool verb;
    float *vel, rho;
    float angle;
    int ix, ih, nh2;
    sf_file velFile, fvx, fvy, snapsf;
    float elapsed, tstart, tend;
    char *antialias;
    const char* mode;

    /* Initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    dip = sf_input ("dip");
    azin = sf_input ("az");
    fvx = sf_input("vx");
    fvy = sf_input("vy");
    velFile = sf_input("vel");

    /* Get dimensions from input */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in inp");
    if (!sf_histint(inp,"n3",&ny)) sf_error("No n3= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"d3",&dy)) sf_error("No d3= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");
    if (!sf_histfloat(inp,"o3",&oy)) sf_error("No o3= in inp");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* Adjoint flag */
    
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    /* Path-integral range */
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001;
    /* Threshold for tail elimination */
    
    if (!sf_getfloat("eps",&eps)) eps = 0.001;
    /* Damper for pi */
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001; 
    /* Damper for t2warp */
    
    /* New axis length */
    if (!sf_getint("pad",&nt2)) nt2=nt;
    /* output time samples */

    //if (!sf_getfloat("vel",&vel)) sf_error("Specify migration velocity");
    /* migration velocity for Kirchhoff */
    
    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    /* Do we use nt from input or the one we specify for the new time axis? */
    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
	
    if (nw < 1 || nw > 3) sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing iline */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing xline */

    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform AzPWD filtering */
    if (!sf_getbool("domod",&domod)) domod=true;
    /* if perform Kirchhoff modeling/migration */
    if (!sf_getbool("dopi",&dopi)) dopi=true;
    /* if perform PI filtering */
 
    if (!sf_getbool("doanisodiff",&doanisodiff)) doanisodiff=true;
    /* if perform anisotropic diffusion regularization */
    if (!sf_getbool("dothr",&dothr)) dothr=true;
    /* if perform sparse regularization */

    if (!sf_getbool("doomp",&doomp)) doomp=false;
    /* OpenMP */

    if (!sf_getbool("snaps",&snaps)) snaps=false;
    /* if do snapshots of outer iterations */

    if (!sf_getint("dsnaps",&dsnaps)) dsnaps=1;
    /* snapshots interval */

    if (!sf_getint("initer",&initer)) initer=2;
    /* inner iterations */
    if (!sf_getint("oniter",&oniter)) oniter=1;
    /* outer iterations */

    if (snaps){
 
	snapsf = sf_output("snapsf"); 

	sf_putint(snapsf,"n1",nt);
	sf_putfloat(snapsf,"d1",dt);
	sf_putfloat(snapsf,"o1",ot);
	sf_putstring(snapsf,"label1","Time");

	sf_putint(snapsf,"n2",nx);
	sf_putfloat(snapsf,"d2",dx);
	sf_putfloat(snapsf,"o2",ox);
	sf_putstring(snapsf,"label2","Distance");

	sf_putint(snapsf,"n3",ny);
	sf_putfloat(snapsf,"d3",dy);
	sf_putfloat(snapsf,"o3",oy);
	sf_putstring(snapsf,"label3","Distance");

	sf_putint(snapsf,"n4",oniter/dsnaps);
	sf_putfloat(snapsf,"d4",dsnaps);
	sf_putfloat(snapsf,"o4",0);
	sf_putstring(snapsf,"label4","Outer Iterations");
	

    } else {snapsf=NULL;}

    if (!sf_getint("niter",&niter)) niter=10;
    /* Anisotropic diffusion: number of conjugate-gradient iterations */
    if (!sf_getint("repeat",&repeat)) repeat=1;
    /* Anisotropic diffusion: number of smoothing iterations */
    if (!sf_getfloat("anisoeps",&anisoeps)) anisoeps=1.;
    /* Anisotropic diffusion: regularization parameter */

    thrflag = sf_getfloat("thr",&thr);
    /* Thresholding level */
    if (!thrflag) sf_error("Need threshold level");
    
    if (thr<0) sf_error("Threshold must be >=0");
    
    mode = sf_getstring("mode");
    /* 'soft', 'hard', 'nng' (default: soft)*/
    if (mode == NULL) mode = "soft";

    n12  = nt*nx*ny;   

    data = sf_floatalloc(n12);
    modl = sf_floatalloc(n12);
    modl0 = sf_floatalloc(n12);

    /* first starting model is zero */
    for(i=0; i<n12; i++){

	modl[i] = 0.0;    	
	modl0[i] = 0.0;

    }

    /* allocate space for dip */
    pp1 = sf_floatalloc(n12);
    pp2 = sf_floatalloc(n12);
    /* allocate space for azimuth */
    az = sf_floatalloc(n12);

    sf_floatread(data,n12,inp);
    
    /* reading iline dip */
    sf_floatread(pp1,n12,dip);
    /* reading xline dip */
    sf_floatread(pp2,n12,dip);
    /* reading azimuth */
    sf_floatread(az,n12,azin);

    /* read structure tensor components */
    vx = sf_floatalloc(n12);
    sf_floatread(vx,n12,fvx);
    vy = sf_floatalloc(n12);
    sf_floatread(vy,n12,fvy);

    /* allocating and reading velocity */
    vel = sf_floatalloc(n12);
    sf_floatread(vel,n12,velFile);

    /*piintface_init(nt,nx,ny, dt,dx,dy, ot,ox,oy, passthr, v_1,v_2,v_3,v_4,
                   eps,epst2, nt2, vel, rho, antialias[0], nw, pp1, pp2,
		   az, nj1, nj2, domod, sm, dopi, doomp, apt, angle); */

    /* outer iterations loop */
    for (j=0; j<oniter; j++){

	if (j == 0) {

		testmig3_init(nt,nx,ny, dt,dx,dy, ot,ox,oy, vel, rho, antialias[0],
				doomp, apt, angle);

		sf_solver(testmig3_lop /* chain of operators loop */,
		      sf_cgstep /* conjugate gradients step */,
		      n12,n12 /* model and data sizes */,
	              modl /* inverted model */,
		      data /* fitted data */,
		      initer /* number of iterations */,
		      //"x0", modl0 /* starting model */,
		      "verb",true,"end");

		sf_cgstep_close();

	} else {

		testmig3_init(nt,nx,ny, dt,dx,dy, ot,ox,oy, vel, rho, antialias[0],
				doomp, apt, angle);

		sf_solver(testmig3_lop /* chain of operators loop */,
		      sf_cgstep /* conjugate gradients step */,
		      n12,n12 /* model and data sizes */,
	              modl /* inverted model */,
		      data /* fitted data */,
		      initer /* number of iterations */,
		      "x0", modl0 /* starting model */,
		      "verb",true,"end");

		sf_cgstep_close();

	}

	for (i3=0; i3<ny; i3++){

		for (i2=0; i2<nx; i2++){
		
			for (i1=0; i1<nt; i1++){

				//if (i1*dt<0.1) modl[i3*nx*nt + i2*nt + i1] = 0.0;	

				if ( ((i1*dt<0.1) && (fabsf(modl[i3*nx*nt + i2*nt + i1]) > 0.01)) && !ch){
					sf_warning("modl[i3*nx*nt + i2*nt + i1]=%f;",modl[i3*nx*nt + i2*nt + i1]);
					ch = true;
				}

			}/* nt */

		}/* nx */
	
	}/* ny */

	/* Writing out snapshot */
	if ((snaps) && (0==j%dsnaps) && (snapsf != NULL)){

		sf_floatwrite(modl,n12,snapsf);

	}

	/* Sparse regularization by thresholding */	
	if (dothr) {
			
    		thr_init(false /* if input is complex */,
             		nx*ny /* number of traces */,
             		nt /* number of time samples */,
	     		thr /* thresholding value */,
	     		mode /*thresholding mode */);

    		thr_lop(modl,modl);
	
	}

	/* Anisotropic diffusion regularization */
	if (doanisodiff){
		
		anisodiffuse_init(nt,nx,ny /* data size */, 
		      vx,vy /* parallel to edges vector components */,
		      niter /* number of iterations */,
		      repeat /* number of smoothing iterations */,
		      anisoeps /* regularization */);

		anisodiffuse_lop(n12,n12,modl,modl0);

	} else {

	/* modl0 = AnisoDIFF( THR(modl) ); Copying to the starting model*/
	for (i=0; i<n12; i++){

		modl0[i] = modl[i];	

		modl[i] = 0.0;

	}}

    }/* outer */

    sf_floatwrite(modl0,n12,out);

    //piintface_close();

    anisodiffuse_close();

    exit(0);

}
