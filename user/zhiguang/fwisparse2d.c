/* Full waveform inversion with sparsity constraint */
/*
  Copyright (C) 2014 University of Texas at Austin
   
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
#include <rsfpwd.h>

#include "fdprep.h"
#include "sparsesolver.h"

bool hermite_false, hermite_true, sparsity;
static float **v, **vnew, **slope, **recloc, *error;
static float **cur_grad, **pre_grad, **cur_dir, **pre_dir;
int n1, n2, n12, npml, pad1, pad2, ns, uts, niter, par;
float d1, d2, pclip;
char *order, *type;
sf_file out;
float *pp, *qq;
int k=0;

void fwisparse_init(sf_file out1, float **v1, float **vnew1, float **slope1, float **recloc1,
		char *order1, int niter1, int uts1, bool hermite_false1, 
		bool hermite_true1, bool sparsity1, float *error1, char *type1, int n11, int n21, int npml1, 
		int pad11, int pad21, int ns1, float d11, float d21, int par1, float pclip1)
/*< initialize >*/
{
	/* output file */
	out=out1;
	/* velocity */
	v=v1;
	vnew=vnew1;
	/* slope */
	slope=slope1;
	recloc=recloc1;

	/* finite difference order */
	order=order1;
	/* number of iteration */
	niter=niter1;
	uts=uts1;

	hermite_false=hermite_false1;
	hermite_true=hermite_true1;
	sparsity=sparsity1;

	/* data misfit array */
	error=error1;
	/* type of seislet transform */
	type=type1;
	/* model dimensions */
	n1=n11;
	n2=n21;
	npml=npml1;
	pad1=pad11;
	pad2=pad21;
	ns=ns1;
	d1=d11;
	d2=d21;
	/* seislet transform accuracy order */
	par=par1;
	pclip=pclip1;

	n12=n1*n2;

	fdprep_order(order);
	if(sparsity) seislet_init(n1,n2,true,true,0.1,par,type[0]);
	if(sparsity) seislet_set(slope);
	if(sparsity){
		pp=sf_floatalloc(n12);
		qq=sf_floatalloc(n12);
	}

	/* storage allocation */
	cur_grad=sf_floatalloc2(n1, n2);
	pre_grad=sf_floatalloc2(n1, n2);
	cur_dir=sf_floatalloc2(n1, n2);
	pre_dir=sf_floatalloc2(n1, n2);

}

void fwisparse_free()
/*< free storage >*/
{
	free(*v); free(v);
	free(*vnew); free(vnew);
	free(*recloc); free(recloc);
	free(error);
	free(*cur_grad); free(cur_grad);
	free(*pre_grad); free(pre_grad);
	free(*cur_dir);  free(cur_dir);
	free(*pre_dir);  free(pre_dir);
	free(pp);
	free(qq);

	if(sparsity) seislet_close();
}

float residual(sf_complex ***obs, sf_complex ***syn, sf_complex ***res)
/*< calculate residual >*/
{
	int i1, i2, is;
	sf_complex temp;
	float misfit=0.0;

	for(is=0; is<ns; is++)
		for(i2=0; i2<n2; i2++)
			for(i1=0; i1<n1; i1++){
				res[is][i2][i1]=sf_cmplx(0.,0.);

				if(recloc[i2][i1]>0.){
					temp=obs[is][i2][i1]-syn[is][i2][i1];
					res[is][i2][i1]=temp;

					misfit+=cabsf(temp)*cabsf(temp);
				}
			}
	return misfit;
}

float gradient(double omega, sf_complex ***f, sf_complex ***obs)
/*< calculate gradient and misfit >*/
{
	int is, i1, i2;
	sf_complex ***syn, ***res;
	float misfit;

	syn=sf_complexalloc3(n1, n2, ns);
	res=sf_complexalloc3(n1, n2, ns);

	/* initialize sparse solver */
	sparse_init(uts, pad1, pad2);

	/* factorize matrix according to different frequencies and models */
	sparse_factor(omega, n1, n2, d1, d2, v, npml, pad1, pad2);

	for(is=0; is<ns; is++)
		for(i2=0; i2<n2; i2++)
			for(i1=0; i1<n1; i1++)
				syn[is][i2][i1]=f[is][i2][i1];

	/* sparse solver for source wavefield */
	sparse_solve(npml, pad1, pad2, syn, hermite_false, ns, uts);

	misfit=residual(obs, syn, res);

	/* sparse solver for receiver wavefield */
	sparse_solve(npml, pad1, pad2, res, hermite_true, ns, uts);

	/* calculate gradient */
	for (i2=0; i2<n2; i2++){
		for (i1=0; i1<n1; i1++){
			cur_grad[i2][i1]=0.;
			for (is=0; is<ns; is++){
				cur_grad[i2][i1] += crealf(conjf(syn[is][i2][i1])*res[is][i2][i1]);
			}
		}
	}

	/* free memory */
	sparse_free(uts);
	free(**syn); free(*syn); free(syn);
	free(**res); free(*res); free(res);

	return misfit;
}

float direction_cg_polak()
/*< calculate conjugate gradient direction >*/
{
	int i1, i2;
	float diff, beta_top, beta_bot, beta;

	beta_top=0.;
	beta_bot=0.;

	for(i2=0; i2<n2; i2++)
		for (i1=0; i1<n1; i1++){
			diff=cur_grad[i2][i1]-pre_grad[i2][i1];
			beta_top += cur_grad[i2][i1]*diff;
			beta_bot += pre_grad[i2][i1]*pre_grad[i2][i1];
		}

	beta=beta_top/beta_bot;

	if (beta < 0.0) beta=0.;

	for (i2=0; i2<n2; i2++){
		for (i1=0; i1<n1; i1++){
			cur_dir[i2][i1]=cur_grad[i2][i1]+beta*pre_dir[i2][i1];
		}
	}
	return beta;
}

void update_model(float alpha)
/*< update model >*/
{
	int i1, i2;
	float dmax=0.;

	for (i2=0; i2<n2; i2++){
		for (i1=0; i1<n1; i1++) {
			if(fabsf(v[i2][i1]) > dmax)
				dmax=fabsf(v[i2][i1]);
		}
	}

	for (i2=0; i2<n2; i2++){
		for (i1=0; i1<n1; i1++){
			vnew[i2][i1]=v[i2][i1]+alpha*cur_dir[i2][i1]/dmax;
		}
	}

}

float forward_operator(double omega, sf_complex ***f, sf_complex ***obs)
/*< forward operator and return data misfit >*/
{
	int is, i2, i1;
	sf_complex ***syn, temp;
	float misfit=0.0;

	syn=sf_complexalloc3(n1, n2, ns);

	/* initialize sparse solver */
	sparse_init(uts, pad1, pad2);

	/* factorize matrix */
	sparse_factor(omega, n1, n2, d1, d2, vnew, npml, pad1, pad2);

	for (is=0; is<ns; is++)
		for (i2=0; i2<n2; i2++)
			for (i1=0; i1<n1; i1++)
				syn[is][i2][i1]=f[is][i2][i1];

	/* sparse solver for source wavefield */
	sparse_solve(npml, pad1, pad2, syn, hermite_false, ns, uts);

	for(is=0; is<ns; is++){
		for (i2=0; i2<n2; i2++){
			for (i1=0; i1<n1; i1++){
				if(recloc[i2][i1]>0.0){
					temp=obs[is][i2][i1]-syn[is][i2][i1];
					misfit += cabsf(temp)*cabsf(temp);
				}
			}
		}
	}

	/* free memory */
	sparse_free(uts);
	
	return misfit;
}

void threshold()
/*< soft thresholding >*/
{
	int m,i;
	float *temp, t, d;

	m=0.5+n1*n2*(1-0.01*pclip);
	if(m<0) m=0;
	if(m>=n12) m=n12-1;

	temp=sf_floatalloc(n12);
	for(i=0; i<n12; i++)
		temp[i]=fabsf(qq[i]);

	t=sf_quantile(m, n12, temp);

	for (i=0; i<n12; i++){
		d=qq[i];
		if(d<-t){
			qq[i]=d+t;
		}else if (d>t){
			qq[i]=d-t;
		}else{
			qq[i]=0.;
		}
	}
}

int fwisparse_exec(double omega, sf_complex ***f, sf_complex ***obs)
/*< nonlinear conjugate gradient iteration >*/
{
	float misfitold=100000000.0, misfit0, misfit1, misfit2, alpha, beta;
	int iter, i2, i1;

	for (iter=1; iter<=niter; iter++){
		sf_warning("Calculating %d out of %d iteration", iter, niter);

		/* calculate gradient */
		misfit0=gradient(omega, f, obs);
		
		/* record data misfit */
		error[k]=misfit0;
		k++;

		if (iter==1){
			for (i2=0; i2<n2; i2++)
				for (i1=0; i1<n1; i1++)
					cur_dir[i2][i1]=cur_grad[i2][i1];
			beta=0.0;
		}else {
			beta=direction_cg_polak();
		}
		sf_warning("Finish direction calculation.");

		/* first test */
		update_model(0.1);
		misfit1=forward_operator(omega, f, obs);
		sf_warning("Finish test1 calculation.");

		/* second test */
		update_model(0.2);
		misfit2=forward_operator(omega, f, obs);
		sf_warning("Finish test2 calculation.");

		/* update model, quadratic fit */
		alpha=(misfit2-4.0*misfit1+3.0*misfit0)/(20.0*(misfit2-2.0*misfit1+misfit0));
		if(alpha<0.0){
			alpha=0.001;
			sf_warning("alpha is smaller than 0.0");
		}

		sf_warning("In iteration %d, alpha = %f, beta=%f, misfit0=%f.", iter, alpha, beta, misfit0);
		sf_warning("Test model alpha=0.1, misfit1=%f, alpha=0.2, misfit2=%f.\n", misfit1, misfit2);

		if (misfit0 > misfitold){
			sf_warning("Terminate at iteration %d.", iter);
			break;
		}else{
			update_model(alpha);

			for(i2=0; i2<n2; i2++){
				for(i1=0; i1<n1; i1++){
					v[i2][i1]=vnew[i2][i1];
					pre_grad[i2][i1]=cur_grad[i2][i1];
					pre_dir[i2][i1]=cur_dir[i2][i1];
				}
			}

			misfitold=misfit0;
		}

		if(sparsity){ /* sparsity regularization */
			/* change from 2D array to 1D */
			for(i2=0; i2<n2; i2++){
				for(i1=0; i1<n1; i1++){
					pp[i2*n1+i1]=v[i2][i1];
				}
			}
			/* adjoint of seislet transform */
			seislet_lop(true, false, n12, n12, qq, pp);
			/* soft thresholding */
			threshold();
			/* seislet transform */
			seislet_lop(false, false, n12, n12, qq, pp);

			/* change from 1D array to 2D */
			for(i2=0; i2<n2; i2++){
				for(i1=0; i1<n1; i1++){
					v[i2][i1]=pp[i2*n1+i1];
				}
			}
		}

		sf_floatwrite(v[0], n1*n2, out);
	} /* end iteration */

	return 0;
}
