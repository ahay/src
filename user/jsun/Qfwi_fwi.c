/* Visco-acoustic FWI */
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
#include "Qfwi_commons.h"
#include "Qfwi_lbfgs.h"
#include "Qfwi_gradient.h"
#include "triutil.h"
/*^*/

void fwi(sf_file Fdat, sf_file Finv, sf_file Fgrad, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_optim optpar, bool verb, int media)
/*< fwi >*/
{
	int iter=0, flag;
	int nz, nx, nzx, nm;
	float fcost;
	float *x, *direction, *grad;
	sf_gradient gradient;
	FILE *fp;

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
	if(fwipar->grad_type==1) {
		if(media==1) gradient=gradient_av;
		else gradient=gradient_v;
		nm=nzx;
		x=array->vv;
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

	optpar->igrad=0;
	optpar->ipair=0;
	optpar->alpha=1.;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
	if(mpipar->cpuid==0){
		l2norm(nm, grad, &optpar->gk_norm);
		print_iteration(fp, iter, optpar);
	}

	/* optimization loop */
	for(iter=0; iter<optpar->niter; iter++){
		if(mpipar->cpuid==0) sf_warning("--------iter=%d---------", iter);

		optpar->ils=0;

		//reverse(nm, grad, direction);
		if(iter==0){
			reverse(nm, grad, direction);
		}else{
			lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
			lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
		} 

		lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
		line_search(nm, x, grad, direction, gradient, optpar, &flag, mpipar->cpuid);
		
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

	if(mpipar->cpuid==0){
		sf_floatwrite(x, nm, Finv);
	}
	if(mpipar->cpuid==0) fclose(fp);
}

void lstri_op(float **dd, float **dwt, float ***ww, float ***mwt, sf_acqui acpar, sf_vec array, sf_pas paspar, bool verb)
/*< ls TRI operator >*/
{
    float **vv;
    int ix,iz,it;

    if (paspar->inv) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
        for         (it=0; it<acpar->nt; it++)
            for     (ix=0; ix<acpar->nx; ix++)
                for (iz=0; iz<acpar->nz; iz++)
                    ww[it][ix][iz] = 0.0f;
        if (NULL!=mwt) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
            for         (it=0; it<acpar->nt; it++)
                for     (ix=0; ix<acpar->nx; ix++)
                    for (iz=0; iz<acpar->nz; iz++)
                        mwt[it][ix][iz] = 1.0f;
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,it)
#endif
        for     (ix=0; ix<acpar->nx; ix++)
            for (it=0; it<acpar->nt; it++)
                dd[ix][it] = 0.0f;
    }

    /* map 1d to 2d */
    vv = (float**) sf_alloc (acpar->nx,sizeof(float*)); 
    vv[0] = array->vv;
    for (ix=1; ix<acpar->nx; ix++) vv[ix] = vv[0]+ix*acpar->nz;

    timerev_init(verb, true, acpar->nt, acpar->nx, acpar->nz, acpar->nb, acpar->rz-acpar->nb, acpar->dt, acpar->dx, acpar->dz, acpar->coef, vv);

    /* calculate model weighting using correlative imaging condition */
    if (paspar->inv && paspar->prec) { 
        if (paspar->ctr) {
            ctimerev(paspar->ngrp,mwt,dd);
            absval(acpar->nz*acpar->nx*acpar->nt,mwt[0][0]);
        } else {
            timerev_lop(true, false, acpar->nz*acpar->nx*acpar->nt, acpar->nt*acpar->nx, mwt[0][0], dd[0]);
            autopow(acpar->nz*acpar->nx*acpar->nt,(float)paspar->ngrp,mwt[0][0]);
        }
        /* smoothing */
        smooth(acpar->nz, acpar->nx, acpar->nt, paspar->rectz, paspar->rectx, paspar->rectt, paspar->repeat, mwt[0][0]);
        /* local normalizaiton */
        swnorm(verb, paspar->sw, acpar->nz, acpar->nx, acpar->nt, paspar->size, paspar->perc, mwt[0][0]);
        /* hard thresholding */
        if (paspar->hard>0) threshold(false, acpar->nz*acpar->nx*acpar->nt, paspar->hard, mwt[0][0]);
    }

    /* apply time-reversal imaging linear operator */
    if (paspar->inv) {
        if (NULL!=dwt) sf_solver(timerev_lop,sf_cgstep,acpar->nz*acpar->nx*acpar->nt,acpar->nt*acpar->nx,ww[0][0],dd[0],paspar->niter,"mwt",mwt[0][0],"wt",dwt[0],"verb",verb,"end");
        else sf_solver(timerev_lop,sf_cgstep,acpar->nz*acpar->nx*acpar->nt,acpar->nt*acpar->nx,ww[0][0],dd[0],paspar->niter,"mwt",mwt[0][0],"verb",verb,"end");
    } else {
        timerev_lop(false, false, acpar->nz*acpar->nx*acpar->nt, acpar->nt*acpar->nx, ww[0][0], dd[0]);
    }
    
    /* close */
    timerev_close();

}

void lstri(sf_file Fdat, sf_file Fmwt, sf_file Fsrc, sf_mpi *mpipar, sf_acqui acpar, sf_vec array, sf_pas paspar, bool verb)
/*< passive source inversion >*/
{
    float **dd, ***ww, ***mwt;

    /*
#ifdef _OPENMP
#pragma omp parallel
    {
        sf_warning("id=%d, nthreads=%d",omp_get_thread_num(),omp_get_num_threads());
    }
#endif */
    dd = sf_floatalloc2(acpar->nt, acpar->nx);
    ww = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
    if (paspar->inv) mwt = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
    else mwt = NULL;
   
    /* read model */
    if (paspar->inv) sf_floatread(dd[0], acpar->nt*acpar->nx, Fdat);
    else sf_floatread(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);

    lstri_op(dd, NULL, ww, mwt, acpar, array, paspar, verb);
    
    if(mpipar->cpuid==0) {
        if (paspar->inv) {
            if (NULL!=Fsrc) sf_floatwrite(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
            if (NULL!=Fmwt) sf_floatwrite(mwt[0][0], acpar->nz*acpar->nx*acpar->nt, Fmwt);
        } else {
            if (NULL!=Fdat) sf_floatwrite(dd[0], acpar->nt*acpar->nx, Fdat);
        }
    }

    /* close */
    free(*dd); free(dd); 
    free(**ww); free(*ww); free(ww);
    if (paspar->inv) { free(**mwt); free(*mwt); free(mwt); }

}

void pfwi(sf_file Fdat, sf_file Finv, sf_file Fgrad, sf_file Fmwt, sf_file Fsrc, sf_mpi *mpipar, sf_sou soupar, sf_acqui acpar, sf_vec array, sf_fwi fwipar, sf_optim optpar, sf_pas paspar, bool verb)
/*< passive fwi >*/
{
	int iter=0, flag;
	int nz, nx, nzx, nm;
	float fcost;
	float *x, *direction, *grad;
	sf_gradient gradient;
	FILE *fp;
        float **dd,**dwt=NULL,***ww,***mwt,***gwt;
        int it,ix,iz;
        /*int wtn1, wtn2, woffn1, woffn2;*/

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
        // JS
        //gradient=gradient_av;
        gradient=gradient_pas_av;
        nm=nzx;
        x=array->vv;

        /* allocate data/source/weight */
        dd = sf_floatalloc2(acpar->nt, acpar->nx);
        ww = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
        gwt= sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
        if (!paspar->onlyvel) {
            mwt = sf_floatalloc3(acpar->nz, acpar->nx, acpar->nt);
            /*
            dwt = sf_floatalloc2(acpar->nt, acpar->nx);

            wtn1=(fwipar->wt1-acpar->t0)/acpar->dt+0.5;
            wtn2=(fwipar->wt2-acpar->t0)/acpar->dt+0.5;
            woffn1=(fwipar->woff1-acpar->r0)/acpar->dr+0.5;
            woffn2=(fwipar->woff2-acpar->r0)/acpar->dr+0.5;
	    residual_weighting(dwt, acpar->nt, acpar->nx, wtn1, wtn2, woffn1, woffn2, !fwipar->oreo);

            sf_file Fdwt;
            Fdwt=sf_output("Fdwt");
            sf_putint(Fdwt,"n1",acpar->nt);
            sf_putint(Fdwt,"n2",acpar->nx);
            sf_floatwrite(dwt[0],acpar->nt*acpar->nx,Fdwt);
            sf_fileclose(Fdwt);
            */
        } else {
            mwt=NULL;
            dwt=NULL;
        }

        /* read data/source */
        // JS
        sf_floatread(dd[0], acpar->nt*acpar->nx, Fdat);
        if (paspar->onlyvel) {
            sf_floatread(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
        } else {
            /* linear inversion of source */
            lstri_op(dd, dwt, ww, mwt, acpar, array, paspar, verb);
        }

        /* calculate gradient mask */
        if (!paspar->onlyvel && paspar->prec && paspar->hidesrc) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
        for         (it=0; it<acpar->nt; it++)
            for     (ix=0; ix<acpar->nx; ix++)
                for (iz=0; iz<acpar->nz; iz++)
                    gwt[it][ix][iz] = mwt[it][ix][iz];
        threshold(true, acpar->nz*acpar->nx*acpar->nt, 0.2, gwt[0][0]);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
        for         (it=0; it<acpar->nt; it++)
            for     (ix=0; ix<acpar->nx; ix++)
                for (iz=0; iz<acpar->nz; iz++)
                    gwt[it][ix][iz] = 1.-gwt[it][ix][iz];
        } else {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(it,ix,iz)
#endif
        for         (it=0; it<acpar->nt; it++)
            for     (ix=0; ix<acpar->nx; ix++)
                for (iz=0; iz<acpar->nz; iz++)
                    gwt[it][ix][iz] = 1.;
        }

	/* initialize */
        // JS
	//gradient_init(Fdat, mpipar, soupar, acpar, array, fwipar, verb);
	gradient_pas_init(dd, ww, gwt, mpipar, acpar, array, fwipar, verb);

	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	/*if(fwipar->onlygrad) return; // program terminates */

        if (!fwipar->onlygrad) {

	if(mpipar->cpuid==0) fp=fopen("iterate.txt","a");
	direction=sf_floatalloc(nm);
	optpar->sk=sf_floatalloc2(nm, optpar->npair);
	optpar->yk=sf_floatalloc2(nm, optpar->npair);

	optpar->igrad=0;
	optpar->ipair=0;
	optpar->alpha=1.;
	optpar->ils=0;
	optpar->fk=fcost;
	optpar->f0=fcost;
        if(mpipar->cpuid==0) {
            l2norm(nm, grad, &optpar->gk_norm);
            print_iteration(fp, iter, optpar);
        }

	/* optimization loop */
        for(iter=0; iter<optpar->niter; iter++){
            if(mpipar->cpuid==0) sf_warning("--------iter=%d---------", iter);

            if (iter%optpar->repeat==0) optpar->alpha=1.;

            optpar->ils=0;

            if(iter==0){
                reverse(nm, grad, direction);
            }else{
                lbfgs_update(nm, x, grad, optpar->sk, optpar->yk, optpar);
                lbfgs_direction(nm, grad, direction, optpar->sk, optpar->yk, optpar);
            } 

            lbfgs_save(nm, x, grad, optpar->sk, optpar->yk, optpar);
            line_search(nm, x, grad, direction, gradient, optpar, &flag, mpipar->cpuid);
            if(mpipar->cpuid==0) {
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

	if(mpipar->cpuid==0){
            sf_floatwrite(x, nm, Finv);
        }

	if(mpipar->cpuid==0) fclose(fp);

        } /* if !onlygrad */

        if (mpipar->cpuid==0 && !paspar->onlyvel) {
            sf_floatwrite(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
            sf_floatwrite(mwt[0][0],acpar->nz*acpar->nx*acpar->nt, Fmwt);
        }

        free(*dd); free(dd); 
        free(**ww); free(*ww); free(ww);
        if(!paspar->onlyvel) { free(**mwt); free(*mwt); free(mwt); }
}
