/* pi operator building wrapping test function flat gaussian weighting smoothing after pi*/
#include <rsf.h>
#include <math.h>
#include "Faddeeva.h"
#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "t2warp.h"
#include "mig2.h"
#include "allp3.h"

//int pi(float * data, int adj);

int main(int argc, char* argv[])
{
    int nt, nt2, nx, i1, i2, n12, i, j;
    bool adj, sm, domod, dd;
    float dt, dt2, dx, ot, ot2, ox, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * data, * output, * datat2, * outputt2, * model;
    sf_file inp, out;
    /* PWD parameters */
    int nw, nj1;
    float *pp, *pwdata;
    sf_file dip;
    /* Kirchhoff params */
    bool hd, verb, doomp, ps, dopi;
    int apt;
    float *v, rho;
    float angle;
    char *antialias;
    int ix;
    sf_file vel;

    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");
    dip = sf_input("dip");

    /* get dimensions from input */
    
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");

    /* get parameters from command line */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    
    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform derivative filtering = PWD */

    if (!sf_getbool("domod",&domod)) domod=true;
    /* if perform modeling via Kirchhoff */

    if (!sf_getbool("doomp",&doomp)) doomp=true;
    /* OMP is forced currently */

    if (!sf_getbool("dopi",&dopi)) dopi=true;
    /* if do pi */

    if (!sf_getbool("ps",&ps)) ps=true;
    /* spherical divergence */

    if (!sf_getbool("dd",&dd)) dd = true;
    /* differentiation in the data domain */
    
    /* Kirchhoff parameters */

    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    /* done in mig2.h */
    //angle = fabsf(tanf(angle*SF_PI/180.0));

    if (!sf_getbool("hd",&hd)) hd = true;
    /* half differentiation */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */
    
    /* path-integral range */
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001;
    /* threshold for tail elimination */
    
    if (!sf_getfloat("eps",&eps)) eps = 0.001;
    /* damper for pi */
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001;
    /* damper for t2warp */
    
    /* new axis length */
    if (!sf_getint("pad",&nt2)) nt2=nt; /* output time samples */
	
    n12 = nt2*nx;   

    data = sf_floatalloc(nt*nx);
    model = sf_floatalloc(nt*nx);
    datat2 = sf_floatalloc(nt2*nx); 
    outputt2 = sf_floatalloc(nt2*nx);
    output = sf_floatalloc(nt*nx);
    
    /* allocate dip */   
    if (sm){

	pp = sf_floatalloc(nt*nx);
        /* allocate space for dip */
	sf_floatread(pp,nt*nx,dip);
	/* read dip */
	pwdata = sf_floatalloc(nt*nx);
	/* allocate space for pwd data */
	
	if (!sf_getint("order",&nw)) nw=1;
        /* [1,2,3] accuracy order */
	
        if (nw < 1 || nw > 3) 
	    sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

        if (!sf_getint("nj1",&nj1)) nj1=1;
        /* antialiasing */

	allpass32d_init(allpass_init(nw, nj1, nt,nx,1, pp));
	/* initialize all-pass-filter structure */

    } else {

	pp = NULL;
	pwdata = NULL;

    }    

    /* allocating and reading velocity */
    v = sf_floatalloc(nt*nx);
    sf_floatread(v,nt*nx,vel);

    if(!adj) {/* fwd */
    
    	if(domod){/* perform modeling */

		/* reading data */
    		sf_floatread(model,nt*nx,inp);

		/* modeling via mig2 */
		mig2_lop(false /* adjoint flag */,
			 false /* addition flag - for x0! do not need it here */,
			 nt,nx /* data size */,
			 dt,dx,
			 ot,ox,
			 data /* zero-offset */,
			 model  /* image */,
			 v /* velocity */,
			 rho /* phase factor */,
			 hd /* half derivative */,
			 antialias[0] /* antialiasing method */,
			 doomp /* perform OpenMP optimization */,
			 apt /* distance aperture in samples */,
			 angle /* angle aperture in degrees */,
			 ps /* spherical divergence */, dd);

    	} else {/* read the data, avert modeling */
    		sf_warning("modelling is disabled");
    		
    		sf_floatread(data,nt*nx,inp);

	}

	if (sm){/* perform PWD */

		allpass32d_lop(adj,false,nt*nx,nt*nx,data,pwdata);
		
		for(i=0;i<nt*nx;i++){
			data[i]=pwdata[i];
		}

	}//PWD flag

    } else {/* adj */
	
    	/* read data */
    	sf_floatread(data,nt*nx,inp);
	
    }/* adj flag */
	
    /* t2warping axis evaluation */ 
    ot2 = ot*ot;
    dt2 = ot+(nt-1)*dt;
    dt2 = (dt2*dt2 - ot2)/(nt2-1);	
    t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);
	
    //sf_warning("t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);\n");
	
    /* compute pi filter */
    //sf_warning("be4 filter");
    flatpifilt_init(nt2, nx, dt2, dx, 0.001, v_1, v_2, v_3, v_4, eps);
    //sf_warning("after filter");

    //sf_warning("pifilt_init(nt2, nx, dt2, dx, v_a, v_b, v_0, beta, eps);\n");

    if(dopi){	

    	if(adj) {

		//sf_warning("be4 the chain: params %d %d %d %d",nt*nx,nt2*nx,nt2*nx,nt*nx);

		sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,adj,false,nt*nx,nt2*nx,nt2*nx,nt*nx,output,data,outputt2,datat2);
	
		//sf_warning("running chain");

    	} else {
	
		sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,data,output,datat2,outputt2);

    	}

    }/* dopi */ else {

	for (i=0;i<nt*nx;i++){

		output[i]=data[i];

	}

    }
	
    if(adj) {

	if (sm){
	
		//sf_warning("performing PWD");
		
		allpass32d_lop(adj,false,nt*nx,nt*nx,pwdata,output);

		for (i=0;i<nt*nx;i++){
				output[i]=pwdata[i];
		}
			
	}

	if (domod) {

		//sf_warning("performing Kirchhoff migration");

		mig2_lop(true /* adjoint flag */,
			 false /* addition flag - for x0! do not need it here */,
			 nt,nx /* data size */,
			 dt,dx,
			 ot,ox,
			 output /* zero-offset */,
			 model  /* image */,
			 v /* velocity */,
			 rho /* phase factor */,
			 hd /* half derivative */,
			 antialias[0] /* antialiasing method */,
			 doomp /* perform OpenMP optimization */,
			 apt /* distance aperture in samples */,
			 angle /* angle aperture in degrees */,
			 ps /* spherical divergence */, dd);

	} else {

		sf_warning("modeling is disabled");
			
		for (i=0;i<nt*nx;i++){
			model[i]=output[i];
		}

	}		

    } /* adj */
	
    //sf_warning("done with output");
	
    if (!adj) {

	sf_floatwrite(output,nt*nx,out);
	
    } else {
	
	sf_floatwrite(model,nt*nx,out);

    }

    exit(0);

}
