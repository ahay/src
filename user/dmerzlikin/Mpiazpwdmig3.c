/* Least-Squares 3D Path-Summation Integral, Azimuthal Plane-Wave Destruction and Kirchhoff Modeling/Migration Chain of Operators*/
#include <rsf.h>
#include <math.h>
#include "piintface.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    int nt, nt2, nx,ny, i1, i2, n12,nt12, i, j;
    bool adj, sm, domod, dopi;
    float dt, dt2, dx,dy, ot, ot2, ox,oy, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * input, * output;
    sf_file inp, out;
    /* PWD parameters */
    int nw, nj1, nj2;
    float *pp1, *pp2, *az;
    sf_file dip, azin;
    /* Kirchhoff params */
    bool half, verb,normalize,debug, doomp;
    int nh, **fold, apt;
    float vel, rho;
    float angle;
    int ix, ih, nh2;
    sf_file velFile;
    float elapsed, tstart, tend;
    char *antialias;

    /* Initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    dip = sf_input ("dip");
    azin = sf_input ("az");

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

    if (!sf_getfloat("vel",&vel)) sf_error("Specify migration velocity");
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

    if (!sf_getbool("doomp",&doomp)) doomp=false;
    /* OpenMP */

    n12  = nt*nx*ny;   

    input  = sf_floatalloc(n12);
    output = sf_floatalloc(n12);

    /* allocate space for dip */
    pp1 = sf_floatalloc(n12);
    pp2 = sf_floatalloc(n12);
    /* allocate space for azimuth */
    az = sf_floatalloc(n12);

    sf_floatread(input,n12,inp);
    
    /* reading iline dip */
    sf_floatread(pp1,n12,dip);
    /* reading xline dip */
    sf_floatread(pp2,n12,dip);
    /* reading azimuth */
    sf_floatread(az,n12,azin);

    piintface_init(nt,nx,ny, dt,dx,dy, ot,ox,oy, passthr, v_1,v_2,v_3,v_4,
                   eps,epst2, nt2, vel, rho, antialias[0], nw, pp1, pp2,
		   az, nj1, nj2, domod, sm, dopi, doomp, apt, angle);

    if (adj == false){

    	piintface_lop(adj,false,n12,n12,input,output);


    } else {

	piintface_lop(adj,false,n12,n12,output,input);

    }

    sf_floatwrite(output,n12,out);

    exit(0);

}

