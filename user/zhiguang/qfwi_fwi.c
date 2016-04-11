/* Acoustic/Visco-acoustic FWI */
/*
 Copyright (C) 2016 The University of Texas at Austin
 
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
#include "qfwi_commons.h"
#include "qfwi_lbfgs.h"
#include "Qfwi_gradient.h"
/*^*/

void fwi(sf_file Fdat, sf_file Finv, sf_file Ferr, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_optim optpar, bool verb, int media)
/*< fwi >*/
{
	int i, iter=0, flag;
	int nz, nx, nzx, nm;
	float fcost, alpha1, alpha2;
	float swap;
	float *x, *direction, *grad, *x1, *x2, *grad1, *grad2, *direction2;
	sf_gradient gradient;
	FILE *fp;
	sf_file Fdir;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
	if(fwipar->grad_type==1) {
		if(media==1) gradient=gradient_av;
		else gradient=gradient_v;
		nm=nzx;
		x=array->vv;
	}else if(fwipar->grad_type==2){
		gradient=gradient_q;
		nm=nzx;
		x=array->tau;
	}else if(fwipar->grad_type==3){
		gradient=gradient_vq;
		nm=2*nzx;
		x=sf_floatalloc(nm);
		for(i=0; i<nzx; i++){
			x[i]=array->vv[i];
			x[nzx+i]=array->tau[i];
		}
		
		if(optpar->opt_type==3) {
			x1=sf_floatalloc(nzx);
			x2=sf_floatalloc(nzx);
			grad1=sf_floatalloc(nzx);
			grad2=sf_floatalloc(nzx);
		}
	}

	/* initialize */
	gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb);

	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	if(fwipar->onlygrad) return; // program terminates 

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");

	direction=sf_floatalloc(nm);
	optpar->sk=sf_floatalloc2(nm, optpar->npair);
	optpar->yk=sf_floatalloc2(nm, optpar->npair);

	optpar->igrad=1;
	optpar->ipair=0;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
	optpar->alpha=1.;
	/* initialize data error vector */
	for(iter=0; iter<optpar->niter+1; iter++){
		optpar->err[iter]=0.;
	}
	take_swap(&swap);
	optpar->err[0]=swap;

	if(mpipar->cpuid==0){
		l2norm(nm, grad, &optpar->gk_norm);
		print_iteration(fp, iter, optpar);
	}

	if(optpar->tangent==1) {
		optpar->niter+=1;
		direction2=sf_floatalloc(nm);
		Fdir=sf_output("Fdir");
		sf_putint(Fdir,"n1",nz);
		sf_putint(Fdir,"n2",nx);
		sf_putint(Fdir,"n3",1);
	}

	/* optimization loop */
	for(iter=0; iter<optpar->niter; iter++){
		if(mpipar->cpuid==0) sf_warning("-------------------iter=%d----------------------", iter+1);
		
		optpar->ils=0;
		if(iter%optpar->repeat==0) optpar->alpha=1.;

		if(optpar->opt_type==2){ // steepest decent
			reverse(nm, grad, direction);
		}else{ //l-bfgs
			if(iter==0){
				reverse(nm, grad, direction);
			}else{
				lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
				lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
				/* tangent scheme */
				if(optpar->tangent==1 && iter==optpar->niter-1){
					/* update sigma1 */
					turnon(1);
					gradient(x, &fcost, grad);
					lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
					for(i=0; i<nm; i++){
						direction2[i]=-direction[i]*optpar->sigma1;
					}
					/* update sigma2 */
					turnon(2);
					gradient(x, &fcost, grad);
					lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
					for(i=0; i<nm; i++){
						direction2[i] -= direction[i]*optpar->sigma2;
					}
					break;
				} // end of tangent scheme
			}
			
			if(optpar->opt_type==3){ // multistep-length gradient
				for(i=0; i<nzx; i++){
					x1[i]=x[i];
					x2[i]=x[i+nzx];
					grad1[i]=grad[i];
					grad2[i]=grad[i+nzx];
				}

				if(mpipar->cpuid==0) sf_warning("-----iter=%d------velocity step length------", iter+1);
				if(iter%optpar->repeat==0) optpar->alpha=1.;
				line_search(nzx, x1, grad1, direction, gradient_v, optpar, &flag, mpipar->cpuid, 31);
				alpha1=optpar->alpha;
				
				update(x);
				if(mpipar->cpuid==0) sf_warning("-----iter=%d------tau step length------", iter+1);
				if(iter%optpar->repeat==0) optpar->alpha=1.;
				line_search(nzx, x2, grad2, direction+nzx, gradient_q, optpar, &flag, mpipar->cpuid, 32);
				alpha2=optpar->alpha;
				
				for(i=0; i<nzx; i++){
					direction[i]=alpha1*direction[i];
					direction[i+nzx]=alpha2*direction[i+nzx];
				}
				
				if(iter%optpar->repeat==0) optpar->alpha=1.;
				if(mpipar->cpuid==0) sf_warning("--iter=%d--velocity and tau step length--", iter+1);
			} // end of multistep-length
			
			lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
		} // end of l-bfgs

		line_search(nm, x, grad, direction, gradient, optpar, &flag, mpipar->cpuid, fwipar->grad_type);
		//if(iter==3) return;
		//
		take_swap(&swap);
		optpar->err[iter+1]=swap;
		
		if(mpipar->cpuid==0){
			l2norm(nm, grad, &optpar->gk_norm);
			print_iteration(fp, iter+1, optpar);
		}

		if(mpipar->cpuid==0 && flag==2){
			fprintf(fp, "Line Search Failed\n");
			break;
		}

		if(mpipar->cpuid==0 && optpar->fk/optpar->f0 < optpar->conv_error){
			fprintf(fp, "Convergence Criterion Reached\n");
			break;
		}
	} // end of iter

	if(mpipar->cpuid==0 && iter==optpar->niter){
		fprintf(fp, "Maximum Iteration Number Reached\n");
	}

	if(mpipar->cpuid==0) sf_floatwrite(x, nm, Finv);
	if(mpipar->cpuid==0) sf_floatwrite(optpar->err, optpar->niter+1, Ferr);
	if(mpipar->cpuid==0 && optpar->tangent==1) sf_floatwrite(direction2, nm, Fdir);
	
	if(mpipar->cpuid==0) fclose(fp);
	if(optpar->tangent==1) free(direction2);

	return;
}
