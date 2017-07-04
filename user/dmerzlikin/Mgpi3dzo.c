/*Gaussian weighting for ZO 3D case*/
#include <rsf.h>
#include <math.h>

#ifdef SF_HAS_COMPLEX_H  
#include "Faddeeva.h"
#endif

//int pi(float * data, int adj);

int main(int argc, char* argv[])
{
    int nt, nx, ny, i1, i2, i3;
    bool verb=true;
    float dt, dx, dy, ot, ox, oy, x, y, t;
    float v_a, v_b, v_0, beta, eps;
    sf_complex * intrace, * outtrace;
    sf_file inp, out;
    double complex alpha, root, u_b, u_a, temp1, temp2, coeff, z; // coefficients for erfi calculation
    

    //MADAGASCAR C API
    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in inp");
    if (!sf_histint(inp,"n3",&ny)) sf_error("No n3= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"d3",&dy)) sf_error("No d3= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");
    if (!sf_histfloat(inp,"o3",&oy)) sf_error("No o3= in inp");
    
    /* get command line parameters */
    if (!sf_getfloat("v_a",&v_a)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_b",&v_b)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_0",&v_0)) beta=0.0; //unweighted pi calculation
    if (!sf_getfloat("beta",&beta)) beta=0.0; 
    if (!sf_getfloat("eps",&eps)) eps = 0.001; // damper for gpi

    dx *= 2.*SF_PI;
    ox *= 2.*SF_PI;

    dy *= 2.*SF_PI;
    oy *= 2.*SF_PI;

    dt *= 2.*SF_PI;
    ot *= 2.*SF_PI;
    
    intrace = sf_complexalloc(nt);
    
    outtrace = sf_complexalloc(nt); 
    
    /* read data currently 3D */
    for (i3=0; i3 < ny; i3++) {
	
        if (verb) sf_warning("wavenumber %d of %d;", i3+1,ny);

	y = oy+i3*dy;
        y *= y;

	for (i2=0; i2 < nx; i2++) {
	    
            x = ox+i2*dx; 
            x *= x;
	    
	    sf_complexread(intrace,nt,inp);
	    
            for (i1=0; i1 < nt; i1++) {
		
            	t = ot+i1*dt; 	

	    	//Path-Integral Analytical Evaluation
	    
	    	//computing coefficients for erfi
	    	alpha = (-1)*((x+eps) + (y+eps))/(16*(t+eps));
			
		root = csqrt(I*alpha - beta);

		//erfi arguments for v_a and v_b 
		u_b = v_b*root + v_0*beta/root;
		u_a = v_a*root + v_0*beta/root;
			
		//integral coefficient	
		coeff = cexp(-beta*v_0*v_0)*cexp(-beta*beta*v_0*v_0/(root*root))/root;
#ifdef SF_HAS_COMPLEX_H  			
		temp1 = coeff*Faddeeva_erfi(u_a,0);
		temp2 = coeff*Faddeeva_erfi(u_b,0);
#else
		temp1 = temp2 = 0;
		sf_error("No C99 complex support");
#endif			
		z = temp2 - temp1;

		z = ((double complex)intrace[i1])*z;

		outtrace[i1] = sf_cmplx(creal(z),cimag(z));

		}// t

		sf_complexwrite (outtrace,nt,out);
	
	}// x

    }// y

    exit(0);
}
