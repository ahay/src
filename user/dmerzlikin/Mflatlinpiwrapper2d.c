/* pi operator building wrapping test function flat gaussian weighting smoothing after pi*/
#include <rsf.h>
#include <math.h>
#include "Faddeeva.h"
#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "t2warp.h"
#include "kirchnew.h"
#include "testmig2.h"

//int pi(float * data, int adj);

int main(int argc, char* argv[])
{
    int nt, nt2, nx, i1, i2, ch, n12, n122, fk;
    bool adj, sm, domod;
    float dt, dt2, dx, ot, ot2, ox, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * data, * datach, * output, * datat2, * outputt2, * smooth, * model;
    sf_file inp, out, pifk;
    /* smoothing variables */
    int nrep, dim, dim1, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM], i0, i, j, nvar;
    bool diff[SF_MAX_DIM], box[SF_MAX_DIM];
    int irep;
    char key[6];
    sf_triangle tr;
    /* kirchhoff params */
    bool half, verb, normalize;
    int nh, **fold, apt;
    float **v, rho, *off;
    float h0, dh, aal, angle;
    char *test;
    int ix, ih, nh2;
    sf_file vel, gather, offset;

    //MADAGASCAR C API
    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");

    /* get parameters from command line */
    if (!sf_getbool("adj",&adj)) adj=false;
    
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

    if (NULL != sf_getstring("gather")) {
	gather = sf_output("gather");
    } else {
	gather = NULL;
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
    
    /* smoothing part - determines various params including dimension along
    which smoothing is performed */
    dim = sf_filedims (inp,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i;
	snprintf(key,6,"diff%d",i+1);
	if (!sf_getbool(key,diff+i)) diff[i]=false;
	/*( diff#=(n,n,...) differentiation on #-th axis )*/
	snprintf(key,5,"box%d",i+1);
	if (!sf_getbool(key,box+i)) box[i]=false;
	/*( box#=(n,n,...) box (rather than triangle) on #-th axis )*/
    }
    
    /* creating parameters for smoothing filter */
    s[0] = 1;
    s[1] = nt;
    //s[2] = ny; for 3D case
    nvar = nt*nx; // 2D 
    // nvar = nt*nx*ny; // 3D
    
    /* to output f-k pi filter */
    if (NULL != sf_getstring("pifk")) {
        pifk = sf_output("pifk");
    } else {
        pifk = NULL;
    }    

    /* if perform derivative filtering */

    if (!sf_getbool("sm",&sm)) sm=true;

    /* if y, do adjoint integration */
    
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001; // threshold for tail elimination
    
    if (!sf_getfloat("eps",&eps)) eps = 0.001; // damper for pi
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001; // damper for t2warp
    
    /* new axis length */
	if (!sf_getint("pad",&nt2)) nt2=nt; /* output time samples */
	
	if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */
	
    //perform modelling
    if (!sf_getbool("domod",&domod)) domod=true;
	
    n12 = nt2*nx;
    
    data = sf_floatalloc(nt*nx);
    model = sf_floatalloc(nt*nx);
    datat2 = sf_floatalloc(nt2*nx); 
    outputt2 = sf_floatalloc(nt2*nx);
    output = sf_floatalloc(nt*nx);

    smooth = NULL;
    
    // allocate space for smoothing
    if(sm){
    	smooth = sf_floatalloc(nt*nx);
    }
    
    // allocating and reading velocity
    v = sf_floatalloc2(nt,nx);
    sf_floatread(v[0],nt*nx,vel);

    if(!adj) {
    
    // read data currently 2D
    	if(domod){
		for (i2=0; i2 < nx; i2++) {
	    		sf_floatread(model+i2*nt,nt,inp);
		}
	mig2_lop(false,half,verb,normalize,nt,nx,nh,fold,apt,data,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);
    	} else
    	{
    		sf_warning("modelling is disabled");
    		for (i2=0; i2 < nx; i2++) {
	    		sf_floatread(data+i2*nt,nt,inp);
		}
	} // internal else
    } else {// adj flag
	
    // read data currently 2D
    for (i2=0; i2 < nx; i2++) {
	    sf_floatread(data+i2*nt,nt,inp);
	}
	
    }// adj flag

    tr = NULL;

	if (adj){	
	if (sm) {
	
		for (j = 0; j < nx; j++) {
			for (i = 0; i < nt; i++) {
		
				smooth[i + j*nt] = data [i + j*nt];
			
			}
		}
	
		// browse through dimensions and smooth
	
		for (i=0; i <= dim1; i++) {
			if (rect[i] <= 1) continue;
			tr = sf_triangle_init (rect[i],n[i],box[i]);
			for (j=0; j < nvar/n[i]; j++) {
				i0 = sf_first_index (i,j,dim1+1,n,s);
				for (irep=0; irep < nrep; irep++) {
					sf_smooth (tr,i0,s[i],diff[i],smooth);	
			    	}
			}		
		}
	
		sf_warning("closing triangle");
		sf_triangle_close(tr);
	
		for (j = 0; j < nx; j++) {
			for (i = 0; i < nt; i++) {
			data[i + j*nt] -= smooth[i + j*nt];	
			}
		}

	}// sm flag
	}
	
	// t2warping axis evaluation 
	ot2 = ot*ot;
	dt2 = ot+(nt-1)*dt;
	dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
	// take in account different output trace length
	t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);
	
	sf_warning("t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);\n");
	
	// compute pi filter
	flatpifilt_init(nt2, nx, dt2, dx, 0.001, v_1, v_2, v_3, v_4, eps);
	
	sf_warning("pifilt_init(nt2, nx, dt2, dx, v_a, v_b, v_0, beta, eps);\n");
	
	if(adj) {

	sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,adj,false,nt*nx,nt2*nx,nt2*nx,nt*nx,output,data,outputt2,datat2);

	} else {
	
	sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,data,output,datat2,outputt2);

	//smoothing	
	if (sm) {
	
		for (j = 0; j < nx; j++) {
			for (i = 0; i < nt; i++) {
		
				smooth[i + j*nt] = output [i + j*nt];
			
			}
		}
	
		// browse through dimensions and smooth

		for (i=0; i <= dim1; i++) {
	    		if (rect[i] <= 1) continue;
	    		tr = sf_triangle_init (rect[i],n[i],box[i]);
	    		for (j=0; j < nvar/n[i]; j++) {
				i0 = sf_first_index (i,j,dim1+1,n,s);
				for (irep=0; irep < nrep; irep++) {
					sf_smooth2 (tr,i0,s[i],diff[i],smooth);		
	    			}		
			}			
		}
	    
		sf_triangle_close(tr);
	
		for (j = 0; j < nx; j++) {
			for (i = 0; i < nt; i++) {
		
				 output [i + j*nt] -= smooth[i + j*nt];
			 
			}
		}

	}// smoothing flag

	}
	
	if(adj) {

		mig2_lop(true,half,verb,normalize,nt,nx,nh,fold,apt,output,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);		

	} // adj flag
	
	sf_warning("done with output");
	
	if (!adj) {
	
		// write
		for (i2=0; i2 < nx; i2++) {
	    	sf_floatwrite(output+i2*nt,nt,out);
	
		}
	} else {
	
		// write
		for (i2=0; i2 < nx; i2++) {
	    	sf_floatwrite(model+i2*nt,nt,out);
	    }
	
	}

    exit(0);
}
