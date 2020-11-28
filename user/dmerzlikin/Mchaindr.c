/* chain diffraction extraction debug version*/
#include <rsf.h>
#include <math.h>
#include "Faddeeva.h"
#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "t2warp.h"
#include "kirchnew.h"
#include "testmig2.h"
#include "allp3.h"

//int pi(float * data, int adj);

int main(int argc, char* argv[])
{
    int nt, nt2, nx, i, j, nx2;
    bool adj, sm, domod;
    float dt, dt2, dx, ot, ot2, ox, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * data, * output, * datat2, * outputt2, * model, * pwdmodel, * outputext, * pwdmodelext;
    sf_file inp, out;
    /* PWD parameters */
    int nw, nj1;
    float *pp = NULL, *pwdata, *pwdataext = NULL;
    sf_file dip,outpwdcheck,outdipcheck;
    /* kirchhoff params */
    bool half, verb,normalize,debug;
    int nh, **fold, apt;
    float **v, rho, *off;
    float h0, dh, aal, angle;
    int ix, ih, nh2;
    sf_file vel, offset;
    /* regularization weight */
    float reg;
    /* files needed for the extended operator */
    float * modelext, * dataext;
    // diffmod - diffraction image from the previous iteration
    // dipim - dip distribution in the image
    sf_file dipim; //diffmod,

    //MADAGASCAR C API
    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");
    dip = sf_input("dip"); //initialize file with dip
    dipim = sf_input("dipim"); //initialize file with image dip
    //diffmod = sf_input("diffmod"); //diffraction model from the previous iteration

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx2)) sf_error("No n2= in inp");
    // half of the traces - actual data
    // half of the traces - model penalization
    nx = (int)(nx2)/2;
    sf_warning("nx=%d nx2=%d \n",nx,nx2);
    sf_warning("nx=%d nx2=%d \n",nx,nx2);
    sf_warning("nx=%d nx2=%d \n",nx,nx2);
    sf_warning("nx=%d nx2=%d \n",nx,nx2);
    sf_warning("nx=%d nx2=%d \n",nx,nx2);
    
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");

    /* get parameters from command line */
    /* adjoint flag */
    if (!sf_getbool("adj",&adj)) adj=false;
    /* if perform derivative filtering = PWD */
    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform modelling via Kirchhoff */
    if (!sf_getbool("domod",&domod)) domod=true;
    /* get regularization parameter */
    if (!sf_getfloat("reg",&reg)) reg=0.0;

    /* debug flag */
    if (!sf_getbool("debug",&debug)){
		
		debug=false;
		outpwdcheck = NULL;
		outdipcheck = NULL;

	} else {

		outpwdcheck = sf_output("outpwd");
    		outdipcheck = sf_output("outdip");

		}
    
    /* kirchhoff parameters */
////////////////////////////////////////////////////////////////////////////////////////

    if (!sf_getbool("normalize",&normalize)) normalize=true;
    /* normalize for the fold */

    if (normalize) {
	fold = sf_intalloc2(nt,nx);
    } else {
	fold = NULL;
    }

    if (adj) {
	if (!sf_histint(inp,"n3",&nh)) sf_error("No n3=");
       
	sf_putint(out,"n3",1);
    } else {
	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
	/* number of offsets (for modeling) */
	
	sf_putint(out,"n3",nh);
    }	

    if (!sf_getfloat("antialias",&aal)) aal = 1.0;
    /* antialiasing */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    angle = fabsf(tanf(angle*SF_PI/180.0));

    if (!sf_getbool("half",&half)) half = true;
    /* if y, the third axis is half-offset instead of full offset */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nh2 = sf_filesize(offset);

	if (nh2 != nh*nx) sf_error("Wrong dimensions in offset");

	off = sf_floatalloc(nh2);	
	sf_floatread (off,nh2,offset);
	sf_fileclose(offset);
    } else {
	if (adj) {
	    if (!sf_histfloat(inp,"o3",&h0)) sf_error("No o3=");
	    if (!sf_histfloat(inp,"d3",&dh)) sf_error("No d3=");
	    sf_putfloat(out,"d3",1.);
	    sf_putfloat(out,"o3",0.);
	} else {
	    if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	    /* offset sampling (for modeling) */
	    if (!sf_getfloat("h0",&h0)) sf_error("Need h0=");
	    /* first offset (for modeling) */
	    sf_putfloat(out,"d3",dh);
	    sf_putfloat(out,"o3",h0);
	}
	
	if (!half) dh *= 0.5;

	off = sf_floatalloc(nh*nx);
	for (ix = 0; ix < nx; ix++) {
	    for (ih = 0; ih < nh; ih++) {
		off[ih*nx+ix] = h0 + ih*dh; 
	    }
	}
	offset = NULL;
    }
////////////////////////////////////////////////////////////////////////////////////////       
    
    /* path-integral range */
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001; // threshold for tail elimination
    
    if (!sf_getfloat("eps",&eps)) eps = 0.001; // damper for pi
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001; // damper for t2warp
    
    /* new axis length */
    if (!sf_getint("pad",&nt2)) nt2=nt; /* output time samples */
	
    /*^^^*/
    /* extended operator arrays */
    dataext = sf_floatalloc(nt*nx);    
    modelext = sf_floatalloc(nt*nx);
    pwdmodel = sf_floatalloc(nt*nx);  
    pwdmodelext = sf_floatalloc(nt*nx);
    outputext = sf_floatalloc(nt*nx);    
    
    /* chain arrays */
    data = sf_floatalloc(nt*nx);
    model = sf_floatalloc(nt*nx);
    datat2 = sf_floatalloc(nt2*nx); 
    outputt2 = sf_floatalloc(nt2*nx);
    output = sf_floatalloc(nt*nx);
    

    pwdata = NULL;
    
    // allocate dip   
    if (sm){

	pp = sf_floatalloc(nt*nx); //allocate space for dip
	sf_floatread(pp,nt*nx,dip); //read dip
	pwdata = sf_floatalloc(nt*nx); //allocate space for pwd data (modelled with kichhoff operator usually)
	pwdataext = sf_floatalloc(nt*nx); //allocate space for pwd data
	
	if (!sf_getint("order",&nw)) nw=1;
        /* [1,2,3] accuracy order */
	
        if (nw < 1 || nw > 3) 
	    sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

        if (!sf_getint("nj1",&nj1)) nj1=1;
        /* antialiasing */

	allpass32d_init(allpass_init(nw, nj1, nt,nx,1, pp)); //initialize all-pass-filter structure

	}    

    // allocating and reading velocity
    v = sf_floatalloc2(nt,nx);
    sf_floatread(v[0],nt*nx,vel);

    if(!adj) {
    
    	if(domod){// perform modelling

                /*^^^*/
		// reading data
                // first we read actual model
    		sf_floatread(model,nt*nx,inp);
                // here we read diffractivity from the
                // previous iteration (should be catenated along
		// the second axis)
                sf_floatread(modelext,nt*nx,inp);

		// modelling via mig2
		//^^^^^
		mig2_lop(false,half,verb,normalize,nt,nx,nh,fold,apt,data,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);
		mig2_lop(false,half,verb,normalize,nt,nx,nh,fold,apt,dataext,v,rho,modelext,off,h0,dh,dx,ot,dt,aal,angle);

    	} else {// just read the data
    		sf_warning("modelling is disabled");
    		
    		sf_floatread(data,nt*nx,inp);
                //sf_floatread(modelext,nt*nx,inp);

	} // internal else

	if (sm){// perform PWD

		//^^^^^
		allpass32d_lop(adj,false,nt*nx,nt*nx,data,pwdata);
                allpass32d_lop(adj,false,nt*nx,nt*nx,dataext,pwdataext);

		if (debug){
			sf_floatwrite(pwdata,nt*nx,outpwdcheck);
			sf_floatwrite(pp,nt*nx,outdipcheck);
		}
		
		//change the address
		for(i=0;i<nt*nx;i++){
			data[i]=pwdata[i];
			dataext[i]=pwdataext[i];
		}

	}//PWD flag

    } else {// adj flag
	
        /*^^^*/
    	// read data currently 2D
    	sf_floatread(data,nt*nx,inp);
        // read extended data - zero PWD over reflection	
        sf_floatread(dataext,nt*nx,inp);
        // since adj PWD applied to zero is zero
	// lets just put zeroes in the correspoding part of the output
	/*for (i1=0;i1<nt*nx;i1++) {
		dataext[i1]=0.0;		
		}*/

    	}// adj flag
	
	// t2warping axis evaluation 
	ot2 = ot*ot;
	dt2 = ot+(nt-1)*dt;
	dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
	// take in account different output trace length
	t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);
	
	sf_warning("t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);\n");
	
	// compute pi filter
	sf_warning("be4 filter");
	flatpifilt_init(nt2, nx, dt2, dx, 0.001, v_1, v_2, v_3, v_4, eps);
	sf_warning("after filter");
	
	sf_warning("pifilt_init(nt2, nx, dt2, dx, v_a, v_b, v_0, beta, eps);\n");
	
	if(adj) {
	sf_warning("be4 the chain: params %d %d %d %d",nt*nx,nt2*nx,nt2*nx,nt*nx);
	sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,adj,false,nt*nx,nt2*nx,nt2*nx,nt*nx,output,data,outputt2,datat2);

	//sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,adj,false,nt*nx,nt2*nx,nt2*nx,nt*nx,outputext,data,outputt2,datat2);

	sf_warning("running chain");

	} else {
	
	sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,data,output,datat2,outputt2);

	sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,dataext,outputext,datat2,outputt2);

	}
	
	if(adj) {

		if (sm){
			sf_warning("performing PWD");

			//^^^^^
			allpass32d_lop(adj,false,nt*nx,nt*nx,pwdata,output);
			//allpass32d_lop(adj,false,nt*nx,nt*nx,pwdataext,outputext);

			//change the address
			for (i=0;i<nt*nx;i++){
				output[i]=pwdata[i];
				//outputext[i]=pwdataext[i];
			}
			
		}

		if (domod) {
			sf_warning("performing Kirchhoff migration");
			//^^^^^
			mig2_lop(true,half,verb,normalize,nt,nx,nh,fold,apt,output,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);
			//mig2_lop(true,half,verb,normalize,nt,nx,nh,fold,apt,outputext,v,rho,modelext,off,h0,dh,dx,ot,dt,aal,angle);
	
		} else {
			sf_warning("changing the address");
			//change the address
			for (i=0;i<nt*nx;i++){
				model[i]=output[i];
			}
		}		

	} // adj flag
	
	sf_warning("done with output");
	
        if (sm){

	// pp = sf_floatalloc(nt*nx); //allocate space for dip
	sf_floatread(pp,nt*nx,dipim); //read dipim

	// we need to re-initialize allpass filter
        // so that it operates in the image domain 
	allpass32d_init(allpass_init(nw, nj1, nt,nx,1, pp)); //initialize all-pass-filter structure

	// all other parameters are supposed to be initialized previously

		if(!adj) {

        		// perform adj pwd in the image domain on the model and diffraction model
        		allpass32d_lop(adj,false,nt*nx,nt*nx,model,pwdmodel);

			for(j=0;j<nx*nt;j++){

				output[j] = output[j] + outputext[j];
				outputext[j] = reg*pwdmodel[j];

			}
	
        	} else {

			allpass32d_lop(adj,false,nt*nx,nt*nx,pwdmodelext,dataext);

			for(j=0;j<nx*nt;j++){
				modelext[j] = model[j];
				model[j] = model[j] + reg*pwdmodelext[j]; 
			}//for

			

		}

	}

        if (!adj) {
	
		// write
	    	sf_floatwrite(output,nt*nx,out);
		// write extended output
		sf_floatwrite(outputext,nt*nx,out);
	
	} else {
	
		// write
		sf_floatwrite(model,nt*nx,out);
		// equivalent to setting reflection component to zero
		sf_floatwrite(modelext,nt*nx,out);

	}

    exit(0);
}
