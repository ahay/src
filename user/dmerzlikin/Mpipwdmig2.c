/* Chain of Path Integral, Plane-Wave Destruction and Kirchhoff migration (based on sfmig2)

works only for zero offset

make sure nh = 1 dh = 1.0 h0 = 0.0 offset file is not used
 
there are flags to disable PWD (Plane-Wave Destruction), P (Path-Integral Filter) and L (Kirchhoff modelling/migration)

no regularization

can be expressed for forward as: data = P PWD L ( reflections + diffractions ) or as a matrix

.                     | reflections  |
| P PWD L   P PWD L | |              | = | data |           
.                     | diffractions |           

can be expressed for adjoint as:

adjoint reflections = L^T PWD^T P^T data

adjoint diffractions = L^T PWD^T P^T data or as a matrix

| reflections  |   | L^T PWD^T P^T |
|              | = |               | | data |
| diffractions |   | L^T PWD^T P^T |

*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#include <math.h>
#include "Faddeeva.h"
#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "t2warp.h"
#include "mig2.h"
#include "allp3.h"

int main(int argc, char* argv[])
{
    int nt, nt2, nx, n12, i, j, nx2, ix;
    bool adj, sm, domod, pi, verb;
    float dt, dt2, dx, ot, ot2, ox, epst2;
    float v_1, v_2, v_3, v_4, eps, passthr;
    float *data, *dataext , *output, *datat2, *outputt2, *model, *pwdmodel, *outputext, *pwdmodelext, *modelext;
    sf_file inp, out, fvel, dip;
    int nw, nj1, apt;
    float *pp, *pwdata, *pwdataext;
    bool hd, ps, doomp, dd;
    float *vel, rho, angle;
    char *antialias;

    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    /* migration/modeling velocity for sfmig2 */
    fvel = sf_input("vel");
    /* dip file is read only if PWD is enabled (sm = true) */

    /* get dimensions from input - depends on the adjoint flag
    search for @@@dimensions@@@ for details */
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&nx2)) sf_error("No n2= in inp");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in inp");
    if (!sf_histfloat(inp,"o1",&ot)) sf_error("No o1= in inp");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in inp");

    /* get parameters from command line */
    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform Plane-Wave destruction (if disabled -> chain = P L) */
    if (!sf_getbool("domod",&domod)) domod=true;
    /* if perform modeling via Kirchhoff (if disabled -> chain = P PWD) */
    if (!sf_getbool("pi",&pi)) pi=true;
    /* if perform Path-Integral filtering (if disabled -> chain = PWD L) */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbose flag */

    if (!sf_getbool("doomp",&doomp)) doomp=true;
    /* OMP flag - currently hard-coded to y */
    
    /* @@@dimensions@@@ - if modeling then we have two models - reflections and diffractions
       [reflections ; diffractions] = [nt*nx ; nt*nx]
       if migration (adjoint=true) then our only input are data - [data]
       [data] = [nt*nx] */

    if (!adj){
        /* [reflections ; diffractions] = [nt*nx ; nt*nx] */
        nx = (int)(nx2)/2;
        if (verb) sf_warning("modeling mode nx=%d nx2=%d nx=nx2/2 \n",nx,nx2);

        sf_putint(out,"n2",nx);

    } else {
        /* [data] = [nt*nx] */
        nx = nx2;
        if (verb) sf_warning("migration mode nx=%d nx2=%d nx=nx2 \n",nx,nx2);

        sf_putint(out,"n2",nx*2);

    }

    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    /* this is done in mig2.h */
    /* angle = fabsf(tanf(angle*SF_PI/180.0)); */

    if (!sf_getbool("hd",&hd)) hd = true;
    /* half derivative */

    if (!sf_getbool("ps",&ps)) ps = true;
    /* amplitude correction */

    if (!sf_getbool("dd",&dd)) dd = true;
    /* differentiation in the data domain */

    //if (!sf_getbool("half",&half)) half = true;
    /* if y, the third axis is half-offset instead of full offset */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* leaky integration constant */       
    
    /* path-integral range:
    path-integral sums images with velocities
    within the range v_2 to v_3;
    v_1 to v_2 and v_3 to v_4 are tapering intervals*/
    if (!sf_getfloat("v_1",&v_1)) sf_error("No integration range specified"); 
    /* no pass velocity */
    if (!sf_getfloat("v_2",&v_2)) sf_error("No integration range specified"); 
    /* first pass velocity */
    if (!sf_getfloat("v_3",&v_3)) sf_error("No integration range specified"); 
    /* second pass velocity */
    if (!sf_getfloat("v_4",&v_4)) sf_error("No integration range specified");  
    /* no pass velocity */
    if (!sf_getfloat("passthr",&passthr)) passthr = 0.001;
    /* threshold for tail elimination */
    if (!sf_getfloat("eps",&eps)) eps = 0.001;
    /* damper for pi */
    if (!sf_getfloat("epst2",&epst2)) epst2 = 0.001;
    /* damper for t2warp */
    
    /* new axis length */
    if (!sf_getint("pad",&nt2)) nt2=nt;
    /* output time samples */
	
    n12 = nt2*nx;   

    /* memory allocation

    suffix "ext" stands for extended model dimensions:

    [reflections ; diffractions] = [model ; modelext].

    */

    if (!adj) {

        dataext = sf_floatalloc(nt*nx);
        outputext = sf_floatalloc(nt*nx);
    
    } else {

        dataext = NULL;
	outputext = NULL;

    }

    modelext = sf_floatalloc(nt*nx);
    pwdmodel = sf_floatalloc(nt*nx);  
    pwdmodelext = sf_floatalloc(nt*nx);
    data = sf_floatalloc(nt*nx);
    model = sf_floatalloc(nt*nx);
    datat2 = sf_floatalloc(nt2*nx); 
    outputt2 = sf_floatalloc(nt2*nx);
    output = sf_floatalloc(nt*nx);
    
    /* dip allocation */   
    if (sm){

	/* dip in the data domain */
    	dip = sf_input("dip");

	pp = sf_floatalloc(nt*nx); // allocate space for dip in data domain
	sf_floatread(pp,nt*nx,dip); // read dip
        /* allocate space for PWD data */
	pwdata = sf_floatalloc(nt*nx);
	pwdataext = sf_floatalloc(nt*nx); // actually not needed for the adjoint mode
	
	if (!sf_getint("order",&nw)) nw=1;
        /* [1,2,3] accuracy order */
	
        if (nw < 1 || nw > 3) 
	    sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

        if (!sf_getint("nj1",&nj1)) nj1=1;
        /* antialiasing */

	/* initialize all-pass-filter structure */
	allpass32d_init(allpass_init(nw, nj1, nt,nx,1, pp));

    } else {

        pwdata = NULL;
	pwdataext = NULL;

    }    

    /* allocating and reading velocity */
    vel = sf_floatalloc(nt*nx);
    sf_floatread(vel,nt*nx,fvel);

    if(!adj) {/* forward */
    
    	if(domod){/* perform modelling */

		/* reading data
                reflections */
    		sf_floatread(model,nt*nx,inp);
                /* diffractions */
                sf_floatread(modelext,nt*nx,inp);
                /* [reflections ; diffractions] */

		/* Kirchhoff modeling (sfmig2) */
		//mig2_lop(false,half,verb,normalize,nt,nx,nh,fold,apt,data,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);
		//mig2_lop(false,half,verb,normalize,nt,nx,nh,fold,apt,dataext,v,rho,modelext,off,h0,dh,dx,ot,dt,aal,angle);
		//mig2_lop (adj,false, nt,nx, dt,dx, ot,ox, trace, out, vel, rho, hd, antialias[0],doomp,apt,angle,ps);

		mig2_lop (false /* adjoint flag */,
			  false /* addition flag - for x0! */,
			  nt, nx /* data size */,
			  dt, dx,
			  ot, ox,
			  data /* zero-offset */,
			  model /* image */,
			  vel /* velocity */,
			  rho /* phase factor */,
			  hd /* half derivative */,
			  antialias[0] /* antialiasing method */,
			  doomp /* perform OpenMP optimization */,
			  apt /* distance aperture in samples */,
			  angle /* angle aperture in degrees */,
			  ps /* spherical divergence */, dd);

		mig2_lop (false /* adjoint flag */,
			  false /* addition flag - for x0! */,
			  nt, nx /* data size */,
			  dt, dx,
			  ot, ox,
			  dataext /* zero-offset */,
			  modelext /* image */,
			  vel /* velocity */,
			  rho /* phase factor */,
			  hd /* half derivative */,
			  antialias[0] /* antialiasing method */,
			  doomp /* perform OpenMP optimization */,
			  apt /* distance aperture in samples */,
			  angle /* angle aperture in degrees */,
			  ps /* spherical divergence */, dd);		

                /* [L reflections ; L diffractions]
                   [data          ; dataext       ] */

		/* if modelling (or other operation is disabled,
		   e.g. domod = false; sm = false; pi = false  )
		   refer to the dimensions and array names at
		   the previous step                            */

    	} else {/* just read the data */
    		if (verb) sf_warning("modelling is disabled");
    		
    		sf_floatread(data,nt*nx,inp);

		sf_floatread(dataext,nt*nx,inp);

	}

	if (sm){/* perform PWD */

		allpass32d_lop(adj,false,nt*nx,nt*nx,data,pwdata);
                allpass32d_lop(adj,false,nt*nx,nt*nx,dataext,pwdataext);
		
		/* change the address */
		for(i=0;i<nt*nx;i++){
			data[i]=pwdata[i];
			dataext[i]=pwdataext[i];
		}

                /* [ PWD L reflections ; PWD L diffractions ]
                   [ data              ; dataext            ] */

	}//PWD flag

    } else {/* adj flag */
	
        /* read data currently 2D */
    	sf_floatread(data,nt*nx,inp);

        /* I will be referring to these data as data_observed */

    }// adj flag

    if (pi) { /* if perform Path-Integral filtering */  
	
	/* t2warping axis evaluation */
    	ot2 = ot*ot;
    	dt2 = ot+(nt-1)*dt;
    	dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
	/* take in account different output trace length */
	t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);
	
	if (verb) sf_warning("t2warp_init(nt,nt2,nx,ot,dt,ot2,dt2,epst2);\n");
	
	/* compute pi filter */
	if (verb) sf_warning("be4 filter");

	flatpifilt_init(nt2, nx, dt2, dx, 0.001, v_1, v_2, v_3, v_4, eps);
	
	if (verb) sf_warning("after filter");
	
    	if (verb) sf_warning("pifilt_init(nt2, nx, dt2, dx, v_a, v_b, v_0, beta, eps);\n");
	
    	if (adj) {
	
        	if (verb) sf_warning("be4 the chain: params %d %d %d %d",nt*nx,nt2*nx,nt2*nx,nt*nx);
		
		sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,adj,false,nt*nx,nt2*nx,nt2*nx,nt*nx,output,data,outputt2,datat2);

        	/* [ output ] <= [ data ]
		[ output ] = [ P^T data_observed ]*/

		if (verb) sf_warning("running chain");

    	} else {
	
        	/* applying Path-Integral to reflections */
		sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,data,output,datat2,outputt2);

        	/* applying Path-Integral to diffractions */
		sf_chain3(t2warp_inv,freqfilt4pi_lop,t2warp,false,false,nt*nx,nt2*nx,nt2*nx,nt*nx,dataext,outputext,datat2,outputt2);

     	   	/* [ P PWD L reflections ; P PWD L diffractions ]
                   [ output              ; outputext            ] 

                   [ output ] <= [ data ]; [ outputext ] <= [ dataext] */

    	}// else
    
    }/* pi flag */ else { /* just changing the addresses */

		if (adj) {
			
			for (i=0;i<nt*nx;i++) {
				/* [ output ] <= [ data ] for adjoint */
				output[i]=data[i];}

		} else {
						
				for (i=0;i<nt*nx;i++) {
				/*[ output ] <= [ data ]; [ outputext ] <= [ dataext] */
				output[i]=data[i];
				outputext[i]=dataext[i];}
				
		}//else

    }/* pi flag disabled */ //else
	
    if(adj) {

        if (sm){

            if (verb) sf_warning("performing PWD");

            allpass32d_lop(adj,false,nt*nx,nt*nx,pwdata,output);

            //change the address
            for (i=0;i<nt*nx;i++){

                output[i]=pwdata[i];

                }

            /* [ output ] = [ PWD^T P^T data_observed ] 
               [ output ] <= [ output ]
            */
			
            }

        if (domod) {
			
            if (verb) sf_warning("performing Kirchhoff migration");

		//mig2_lop(true,half,verb,normalize,nt,nx,nh,fold,apt,output,v,rho,model,off,h0,dh,dx,ot,dt,aal,angle);

		mig2_lop (true /* adjoint flag */,
			  false /* addition flag - for x0! */,
			  nt, nx /* data size */,
			  dt, dx,
			  ot, ox,
			  output /* zero-offset */,
			  model /* image */,
			  vel /* velocity */,
			  rho /* phase factor */,
			  hd /* half derivative */,
			  antialias[0] /* antialiasing method */,
			  doomp /* perform OpenMP optimization */,
			  apt /* distance aperture in samples */,
			  angle /* angle aperture in degrees */,
			  ps /* spherical divergence */, dd);

            /* [ model ] = [ L^T PWD^T P^T data_observed ] 

               [ model ] <= [ output ]             */
	
        } else {
         
            if (verb) sf_warning("changing the address");

            for (i=0;i<nt*nx;i++){

                model[i]=output[i];

            }

        }		

    } // adj flag
	
    if (verb) sf_warning("done with output");

    if (adj) {

        for(j=0;j<nx*nt;j++){
	    
            modelext[j] = model[j];

            /* [ model (reflections) ] = [ L^T PWD^T P^T data_observed ]
               [ modelext (diffractions) ] = [ L^T PWD^T P^T data_observed ]
            */ 

        }

    }

    else { /* not adjoint */

        for(j=0;j<nx*nt;j++){

            output[j] = output[j] + outputext[j];
            
            /* [ P PWD L reflections ; P PWD L diffractions ]
               [ output              ; outputext            ] 

            [ output ] = [ P PWD L reflections + P PWD L diffractions ] */

        }

    }

    if (!adj) {
	
        sf_floatwrite(output,nt*nx,out);
	
    } else {
	

        sf_floatwrite(model,nt*nx,out);
	
        sf_floatwrite(modelext,nt*nx,out);

    }

    exit(0);
}
