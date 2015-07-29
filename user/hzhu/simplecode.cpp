#include <iostream>
#include "rsf.hh" 

int 
main( int argc, char* argv[]) 
{

	// read in observations 


	// read in starting model


	// read in frequency group and iteration number 


	// set up acquisition 
	gen_acquisition(srcx,srcy,srcz,recx,recy,recz);

	// loop over frequency 
	for (int iw=0; iw<nw; iw++ ) { 
		freq=ow+iw*dw; // frequency 

		// loop over cg iterations
		for (int iter=0; iter< niter; iter++ ) {

			// helm solver forward 
			helm_solver(model,freq,solution,adjoint);

			// generate adjoint sources and misfits 
			gen_adjsrc(solution,obs,adjreal,adjimag,misfit,recx,recy,recz,nx,ny,nz,nsrc,nrec);

			// helm solver adjoint 
			helm_solver(model,freq,adjreal,adjimag,kernel,adjoint);

			// smooth kernel 
			smooth_kernel(kernel,kernel_smooth,nx,ny,nz,sigma2); 

			curker=kernel_smooth;

			// calculate direction 
			if ( iter ==0 ) {
				calc_gradient_sd(curker,curdir,nx,ny,nz);
			else 
				calc_gradient_cg(curker,preker,predir,curdir,nx,ny,nz); 
			}

			// line search
			misfit0=misfit;
			for (istep=0, step=1.0; istep < nstep; istep++, step*=0.5) {
			
				update_model(model,curdir,model_test,step,nx,ny,nz); 
				
				// helm solver forward
				helm_solver(model_test,freq,solution,adjoint);

				// calc misfit value 
				gen_adjsrc(solution,obs,adjreal,adjimag,misfit1,recx,recy,recz,nx,ny,nz,nsrc,nrec);
				if (misfit1>misfit0) { 
					alpha=step;	
					break;
				}
				misfit0=misfit1;
			}	

			// model update 
			update_model(model,curdir,model_new,alpha,nx,ny,nz);
			preker=curker;
			predir=curdir;
			model=model_new;

		}
	}

	return 0;
}


