/* pi operator building wrapping test function flat gaussian weighting smoothing after pi*/
#include <rsf.h>
#include <math.h>
#include "Faddeeva.h"
#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "t2warp.h"
#include "kirchnew.h"

//int pi(float * data, int adj);

int main(int argc, char* argv[])
{
    int nt, nt2, nx, i1, i2, n122=0;
    bool adj, sm, domod;
    float dt, dt2, dx, ot, ot2, ox, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float * data, * output, * datat2, * outputt2, * smooth=NULL, * model;
    sf_file inp, out;
    /* smoothing variables */
    int nrep, dim, dim1, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM], i0, i, j, nvar;
    bool diff[SF_MAX_DIM], box[SF_MAX_DIM];
    int irep;
    char key[6];
    sf_triangle tr=NULL;
    /* kirchhoff params */
    bool hd;
    int sw;
    float *vrms, v0;
    char *test;
    sf_file vel;

    //MADAGASCAR C API
    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");
    
    /* kirchhoff parameters */
    if (!sf_getbool("hd",&hd)) hd=true;
    if (!sf_getbool("domod",&domod)) domod=true;
    /* if y, apply half-derivative filter */
    if (!sf_getint("sw",&sw)) sw=0;
    /* if > 0, select a branch of the antialiasing operation */
    
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
    
    /* get parameters from command line */
    if (!sf_getbool("adj",&adj)) adj=false;

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
	
    data = sf_floatalloc(nt*nx);
    model = sf_floatalloc(nt*nx);
    datat2 = sf_floatalloc(nt2*nx); 
    outputt2 = sf_floatalloc(nt2*nx);
    output = sf_floatalloc(nt*nx);

    /* allocate space for smoothing */
    if(sm){
    	smooth = sf_floatalloc(nt*nx);
    }
    
    /* allocating and reading velocity */
    vrms = sf_floatalloc(nt);
    if (NULL != (test = sf_getstring("velocity"))) { 
	/* velocity file */
	free(test);
	vel = sf_input("velocity");
	sf_floatread(vrms,nt,vel);
	sf_fileclose(vel);
    } else {
	if (!sf_getfloat("v0",&v0)) sf_error("Need velocity= or v0=");
	/* constant velocity (if no velocity=) */
	for (i1=0; i1 < nt; i1++) {
	    vrms[i1] = v0;
	}
    }
    
    if(domod){// if perform modelling
	/* initialize kirchhoff */
    	kirchnew_init (vrms, ot, dt, dx, nt, nx, sw, true, hd);
    	n122 = nt*nx;
    }

    if(!adj) {
    
    /* read data currently 2D */
    	if(domod){
		for (i2=0; i2 < nx; i2++) {
	    		sf_floatread(model+i2*nt,nt,inp);
		}
	kirchnew_lop (false,false,n122,n122,model,data);
    	} else
    	{
    		sf_warning("modelling is disabled");
    		for (i2=0; i2 < nx; i2++) {
	    		sf_floatread(data+i2*nt,nt,inp);
		}
	} // internal else
    } else {// adj flag
	
    /* read data currently 2D */
    for (i2=0; i2 < nx; i2++) {
	    sf_floatread(data+i2*nt,nt,inp);
	}
	
	}

	if (adj){	
	if (sm) {
	
		for (j = 0; j < nx; j++) {
			for (i = 0; i < nt; i++) {
		
				smooth[i + j*nt] = data [i + j*nt];
			
			}
		}
	
		/* browse through dimensions and smooth*/
	
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
	
	/* t2warping axis evaluation */
	ot2 = ot*ot;
	dt2 = ot+(nt-1)*dt;
	dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
	/* take in account different output trace length */
	t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);
	
	sf_warning("t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);\n");
	
	/* compute pi filter */
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
	
		/* browse through dimensions and smooth*/

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

	kirchnew_lop (true,false,n122,n122,model,output);
	sf_warning("chain is done;\n");
	} // adj flag
	
	sf_warning("done with output");
	
	if (!adj) {
	
		/* write */
		for (i2=0; i2 < nx; i2++) {
	    	sf_floatwrite(output+i2*nt,nt,out);
	
		}
	} else {
	
		/* write */
		for (i2=0; i2 < nx; i2++) {
	    	sf_floatwrite(model+i2*nt,nt,out);
	    }
	
	}

    exit(0);
}
