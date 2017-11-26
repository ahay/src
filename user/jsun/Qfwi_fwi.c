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

void lstri(sf_file Fdat, sf_file Fmwt, sf_file Fsrc, sf_mpi *mpipar, sf_acqui acpar, sf_vec array, sf_pas paspar, bool verb)
/*< passive source inversion >*/
{
    float **dd, ***ww, ***mwt;
    int nturn, iturn, is, rdn;
    char filename[20]="tempbin",srdn[10];
    FILE *temp;
    MPI_Comm comm=MPI_COMM_WORLD;

    if (mpipar->cpuid==0) {
        srand(time(NULL));
        rdn = rand()%1000000000;
        sprintf(srdn,"%d",rdn);
        strcat(filename,srdn);
    }
    MPI_Bcast(filename, 20, MPI_CHAR, 0, comm);
    if(verb && mpipar->cpuid==0) sf_warning("filename=%s",filename);

    temp=fopen(filename, "wb+");

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

    if(acpar->ns%mpipar->numprocs==0) nturn=acpar->ns/mpipar->numprocs;
    else nturn=acpar->ns/mpipar->numprocs+1;

    for(iturn=0; iturn<nturn; iturn++){
        is=iturn*mpipar->numprocs+mpipar->cpuid;
        if(is<acpar->ns){
            if (paspar->inv) {
                /* read data */
                sf_seek(Fdat, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                sf_floatread(dd[0], acpar->nt*acpar->nx, Fdat);
            } else {
                /* read source */
                sf_seek(Fsrc, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                sf_floatread(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
            }

            /* do the computation */
            lstri_op(dd, NULL, ww, mwt, acpar, array, paspar, verb);

            if (paspar->inv) {
                /* write source */
                fseeko(temp, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                fwrite(ww[0][0], sizeof(float), acpar->nz*acpar->nx*acpar->nt, temp);
                if (NULL!=Fmwt && is==0) sf_floatwrite(mwt[0][0], acpar->nz*acpar->nx*acpar->nt, Fmwt);
            } else {
                /* write data */
                fseeko(temp, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                fwrite(dd[0], sizeof(float), acpar->nt*acpar->nx, temp);
            }

        } /* if is<ns */
    }
    fclose(temp);
    MPI_Barrier(comm);
    
    if(mpipar->cpuid==0) {
        temp=fopen(filename, "rb");
        if (paspar->inv) {
            for(is=0; is<acpar->ns; is++){
                fseeko(temp, is*acpar->nz*acpar->nx*acpar->nt*sizeof(float), SEEK_SET);
                fread(ww[0][0], sizeof(float), acpar->nz*acpar->nx*acpar->nt, temp);
                sf_floatwrite(ww[0][0], acpar->nz*acpar->nx*acpar->nt, Fsrc);
            }
        } else {
            for(is=0; is<acpar->ns; is++){
                fseeko(temp, is*acpar->nt*acpar->nx*sizeof(float), SEEK_SET);
                fread(dd[0], sizeof(float), acpar->nt*acpar->nx, temp);
                sf_floatwrite(dd[0], acpar->nt*acpar->nx, Fdat);
            }
        }
        fclose(temp);
        remove(filename);
    }
    MPI_Barrier(comm);

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

	nz=acpar->nz;
	nx=acpar->nx;
	nzx=nz*nx;

	/* gradient type */
        gradient=gradient_pas_av;
        nm=nzx;
        x=array->vv;

	/* initialize */
	gradient_pas_init(Fdat, Fsrc, Fmwt, mpipar, soupar, acpar, array, fwipar, paspar, verb);

	/* calculate first gradient */
	grad=sf_floatalloc(nm);
	gradient(x, &fcost, grad);

	/* output first gradient */
	if(mpipar->cpuid==0) sf_floatwrite(grad, nm, Fgrad);

	if(fwipar->onlygrad) return; /* program terminates */

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

}
