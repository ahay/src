/* L-BFGS optimization and line search based on Wolfe condition */
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
/*^*/

void lbfgs_save(int n, float *x, float *grad, float **sk, float **yk, sf_optim opt)
/*< save current model and gradient >*/
{
	int i;
	if(opt->ipair < opt->npair){
		copy(n, x, sk[opt->ipair]);
		copy(n, grad, yk[opt->ipair]);
		opt->ipair += 1;
	}else{
		for(i=0; i<opt->npair-1; i++){
			copy(n, sk[i+1], sk[i]);
			copy(n, yk[i+1], yk[i]);
		}
		copy(n, x, sk[opt->npair-1]);
		copy(n, grad, yk[opt->npair-1]);
	}
}

void lbfgs_update(int n, float *x, float *grad, float **sk, float **yk, sf_optim opt)
/*< update current sk and yk >*/
{
	int i, j;
	j=opt->ipair-1;
	for(i=0; i<n; i++){
		sk[j][i]=x[i]-sk[j][i];
		yk[j][i]=grad[i]-yk[j][i];
	}
}

void lbfgs_direction(int n, float *grad, float *r, float **sk, float **yk, sf_optim opt)
/*< calculate search direction >*/
{
	int i, j;
	float *rho, *q, *alpha, tmp, tmp1, gamma, beta;

	// safeguard
	l2norm(n, sk[opt->ipair-1], &tmp);
	l2norm(n, yk[opt->ipair-1], &tmp1);
	if(tmp==0. || tmp1==0.){
		reverse(n, grad, r);
		return;
	}
	
	q=sf_floatalloc(n);
	rho=sf_floatalloc(opt->ipair);
	alpha=sf_floatalloc(opt->ipair);

	copy(n, grad, q);
	
	// first loop
	for(i=opt->ipair-1; i>=0; i--){
		
		// calculate rho
		dot_product(n, yk[i], sk[i], &tmp);  
		rho[i]=1./tmp;

		dot_product(n, sk[i], q, &tmp);
		alpha[i]=rho[i]*tmp;
		for(j=0; j<n; j++)
			q[j] -= alpha[i]*yk[i][j];
	}

	// initial Hessian
	dot_product(n, yk[opt->ipair-1], yk[opt->ipair-1], &tmp);
	gamma=1./tmp/rho[opt->ipair-1];
	for (j=0; j<n; j++){
		r[j]=gamma*q[j];
	}

	// second loop
	for(i=0; i<opt->ipair; i++){
		dot_product(n, yk[i], r, &tmp);
		beta=tmp*rho[i];
		tmp=alpha[i]-beta;
		for(j=0; j<n; j++)
			r[j] += tmp*sk[i][j];
	}

	// opposite direction of H^(-1)*grad(f)
	for(j=0; j<n; j++)
		r[j]=-r[j];

	// deallocate variables
	free(q);
	free(alpha);
	free(rho);
}

void clip(float *x, int n, float min, float max)
/*< clip data >*/
{
	int i;

	for(i=0; i<n; i++){
		if(x[i]<min)
			x[i]=min;
		if(x[i]>max)
			x[i]=max;
	}
}

void line_search(int n, float *x, float *grad, float *direction, sf_gradient gradient, sf_optim opt, int *flag, int cpuid, int grad_type)
/*< line search (Wolfe condition) >*/
{
	int i, j;
	float m1, m2, m3, fcost, alpha1=0., alpha2=0.;
	float *xk;

	xk=sf_floatalloc(n);
	copy(n, x, xk);
	dot_product(n, grad, direction, &m1);
	m2=m1*opt->c2;
	m1=m1*opt->c1;
	
	for(i=0; i<opt->nls; i++){
		
		opt->ils += 1;
		for(j=0; j<n; j++)
			x[j] =xk[j] + opt->alpha*direction[j];

		/* model constraints */
		if(grad_type==1 || grad_type==31){
			clip(x, n, opt->v1, opt->v2);
		}else if(grad_type==2 || grad_type==32){
			clip(x, n, opt->tau1, opt->tau2);
		}else if(grad_type==3){
			clip(x, n/2, opt->v1, opt->v2);
			clip(x+n/2, n/2, opt->tau1, opt->tau2);
		}

		gradient(x, &fcost, grad);
		opt->igrad += 1;
		dot_product(n, grad, direction, &m3);
		
		if(cpuid==0){
			sf_warning("line search i=%d",i+1);
			sf_warning("grad_type=%d opt->v2=%f opt->tau2=%f", grad_type, opt->v2, opt->tau2);
			sf_warning("alpha1=%f alpha2=%f alpha=%f",alpha1, alpha2, opt->alpha);
			sf_warning("fcost=%f fk=%f fk+c1*alpha*m1=%f m3=%f c2*m1=%f ",fcost, opt->fk, opt->fk+opt->alpha*m1, m3, m2);
		}

		//if(counter1==4 && i==0) return;

		if(opt->opt_type==3 && grad_type ==32){ // Armijo condition
			if(fcost <= opt->fk + opt->alpha*m1){
				if(grad_type!=31 && grad_type!=32) opt->fk=fcost;
				*flag=0;
				break;
			}else{
				alpha2=opt->alpha;
				opt->alpha=0.5*(alpha1+alpha2);
			}
		}else{ // Wolfe condition 
			
			if(fcost <= opt->fk + opt->alpha*m1 && m3 >= m2){
				if(grad_type!=31 && grad_type!=32) opt->fk=fcost;
				//	opt->fk=fcost;
				*flag=0;
				break;
			}else if (fcost > opt->fk + opt->alpha*m1){
				alpha2=opt->alpha;
				opt->alpha=0.5*(alpha1+alpha2);
			}else{
				alpha1=opt->alpha;
				if(alpha2 == 0.)
					opt->alpha *= opt->factor;
				else
					opt->alpha = 0.5*(alpha1+alpha2);
			}
			if(cpuid==0){
				sf_warning("alpha1=%f alpha2=%f alpha=%f",alpha1, alpha2, opt->alpha);
			}
		} // end of Wolfe condition

	} // end of line search
	
	if(i==opt->nls){
		if(fcost <= opt->fk)
			*flag=1;
		else
			*flag=2;
	}

	free(xk);
}
