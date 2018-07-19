/* 3D Path-Summation Integral Operator as a Linear Filter*/
#include <rsf.h>
#include <math.h>
#include "tdpi.h"
#include "Faddeeva.h"
#include "fft3.h"
#include "t2warp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    int nt, nt2, nx,ny, i1, i2, n12,nt12, i, j;
    bool adj, sm, domod;
    float dt, dt2, dx,dy, ot, ot2, ox,oy, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * data, * output, * datat2, * outputt2, * model;
    sf_file inp, out;
    /* PWD parameters */
    int nw, nj1;
    float *pp, *pwdata;
    sf_file dip,outpwdcheck,outdipcheck;
    /* kirchhoff params */
    bool half, verb,normalize,debug;
    int nh, **fold, apt;
    float **v, rho, *off;
    float h0, dh, aal, angle;
    int ix, ih, nh2;
    sf_file vel, gather, offset;
    float elapsed, tstart, tend;

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

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    /* path-integral range */
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001; // threshold for tail elimination
    
    if (!sf_getfloat("eps",&eps)) eps = 0.001; // damper for pi
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001; // damper for t2warp
    
    /* new axis length */
    if (!sf_getint("pad",&nt2)) nt2=nt; /* output time samples */

    n12  = nt*nx*ny; 	
    nt12 = nt2*nx*ny;   

    data = sf_floatalloc(n12);
    model = sf_floatalloc(n12);
    datat2 = sf_floatalloc(nt12); 
    outputt2 = sf_floatalloc(nt12);
    output = sf_floatalloc(n12);

    if (adj == false){

    	sf_floatread(model,n12,inp);

	/* I will add Kirchhoff modeling */
	for(i=0; i<n12; i++){

		data[i] = model[i];

	}

    }

    if (adj == true){
	
    	sf_floatread(data,n12,inp);
	
    }/* adj flag */
	
    /* t2warping axis evaluation */
    ot2 = ot*ot;
    dt2 = ot+(nt-1)*dt;
    dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
    /* take in account different output trace length */
    t2warp_init(nt,nt2,nx*ny,ot,dt,ot2,dt2,epst2);

    #ifdef _OPENMP
    tstart = omp_get_wtime();
    #endif
	
    /* compute pi filter */
    tdpi_init(nt2, nx,ny, dt2, dx,dy, 0.001, v_1, v_2, v_3, v_4, eps);
	
    if(adj == true) {

	sf_chain3(t2warp_inv,tdpi_lop,t2warp,adj,false,nt*nx*ny,nt2*nx*ny,nt2*nx*ny,nt*nx*ny,output,data,outputt2,datat2);

    } else {
	
	sf_chain3(t2warp_inv,tdpi_lop,t2warp,false,false,nt*nx*ny,nt2*nx*ny,nt2*nx*ny,nt*nx*ny,data,output,datat2,outputt2);

    }
    
    #ifdef _OPENMP
    tend = omp_get_wtime();
    elapsed = tend - tstart;
    sf_warning("elapsed=%f [s]",elapsed);
    #endif

    if (adj == true ) {

	/* I will add Kirchhoff modeling */
	for(i=0; i<n12; i++){

		model[i] = output[i];

	}
    	

    }

    if (adj == false) {
	
	sf_floatwrite(output,nt*nx*ny,out);
	
    } else {
	
	sf_floatwrite(model,nt*nx*ny,out);

    }

    exit(0);

}
