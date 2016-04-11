/* Acoustic/Visco-acoustic Forward Modeling, FWI, and RTM */
/*
 Copyright (C) 2016 University of Texas at Austin
 
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
#include <mpi.h>
#include "Qfwi_modeling.h"
#include "qfwi_fwi.h"
#include "Qfwi_rtm.h"

int main(int argc, char* argv[])
{
	int function, media;
	bool verb;

	sf_mpi mpipar;
	sf_sou soupar;
	sf_acqui acpar;
	sf_vec array;
	sf_fwi fwipar;
	sf_optim optpar=NULL;

	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Fq, Ftau, Fw, Fdat, Fimg, Finv=NULL, Ferr=NULL, Fgrad;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &mpipar.cpuid);
	MPI_Comm_size(comm, &mpipar.numprocs);

	sf_init(argc, argv);

	Fv=sf_input("Fvel"); /* velocity model */
	Fq=sf_input("Fq"); /* quality factor */
	Ftau=sf_input("Ftau"); /* tau which determines the magnitude of Q */
	Fw=sf_input("Fwavelet"); /* wavelet */

	soupar=(sf_sou)sf_alloc(1, sizeof(*soupar));
	acpar=(sf_acqui)sf_alloc(1, sizeof(*acpar));
	array=(sf_vec)sf_alloc(1, sizeof(*array));

	/* parameters I/O */
	if(!sf_getint("media", &media)) media=2;
	/* if 1, acoustic media; if 2, visco-acoustic media */
	if(!sf_getint("function", &function)) function=2;
	/* if 1, forward modeling; if 2, FWI; if 3, RTM */

	if(!sf_histint(Fv, "n1", &acpar->nz)) sf_error("No n1= in Fv");
	if(!sf_histint(Fv, "n2", &acpar->nx)) sf_error("No n2= in Fv");
	if(!sf_histfloat(Fv, "d1", &acpar->dz)) sf_error("No d1= in Fv");
	if(!sf_histfloat(Fv, "d2", &acpar->dx)) sf_error("No d2= in Fv");
	if(!sf_histfloat(Fv, "o1", &acpar->z0)) sf_error("No o1= in Fv");
	if(!sf_histfloat(Fv, "o2", &acpar->x0)) sf_error("No o2= in Fv");
	if(!sf_histint(Fw, "n1", &acpar->nt)) sf_error("No n1= in Fw");
	if(!sf_histfloat(Fw, "d1", &acpar->dt)) sf_error("No d1= in Fw");
	if(!sf_histfloat(Fw, "o1", &acpar->t0)) sf_error("No o1= in Fw");

	if(!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
	if(!sf_getint("nb", &acpar->nb)) acpar->nb=100; /* boundary width */
	if(!sf_getfloat("coef", &acpar->coef)) acpar->coef=0.002; /* absorbing boundary coefficient */

	if(!sf_getint("acqui_type", &acpar->acqui_type)) acpar->acqui_type=1;
	/* if 1, fixed acquisition; if 2, marine acquisition; if 3, symmetric acquisition */
	if(!sf_getint("ns", &acpar->ns)) sf_error("shot number required"); /* shot number */
	if(!sf_getfloat("ds", &acpar->ds)) sf_error("shot interval required"); /* shot interval */
	if(!sf_getfloat("s0", &acpar->s0)) sf_error("shot origin required"); /* shot origin */
	if(!sf_getint("sz", &acpar->sz)) acpar->sz=5; /* source depth */
	if(!sf_getint("nr", &acpar->nr)) acpar->nr=acpar->nx; /* number of receiver */
	if(!sf_getfloat("dr", &acpar->dr)) acpar->dr=acpar->dx; /* receiver interval */
	if(!sf_getfloat("r0", &acpar->r0)) acpar->r0=acpar->x0; /* receiver origin */
	if(!sf_getint("rz", &acpar->rz)) acpar->rz=5; /* receiver depth */

	if(!sf_getfloat("f0", &acpar->f0)) sf_error("reference frequency required"); /* reference frequency */
	if(!sf_getint("interval", &acpar->interval)) acpar->interval=1; /* wavefield storing interval */

	if(!sf_getfloat("fhi", &soupar->fhi)) soupar->fhi=0.5/acpar->dt; /* high frequency in band, default is Nyquist */
	if(!sf_getfloat("flo", &soupar->flo)) soupar->flo=0.; /* low frequency in band, default is zero */
	if(!sf_getint("rectx", &soupar->rectx)) soupar->rectx=2; /* source smoothing in x */
	if(!sf_getint("rectz", &soupar->rectz)) soupar->rectz=2; /* source smoothing in z */

	/* get prepared */
	preparation(Fv, Fq, Ftau, Fw, acpar, soupar, array);

	if(function == 1){ // forward modeling
		Fdat=sf_output("output"); /* shot data */

		/* dimension set up */
		sf_putint(Fdat, "n1", acpar->nt);
		sf_putfloat(Fdat, "d1", acpar->dt);
		sf_putfloat(Fdat, "o1", acpar->t0);
		sf_putstring(Fdat, "label1", "Time");
		sf_putstring(Fdat, "unit1", "s");
		sf_putint(Fdat, "n2", acpar->nr);
		sf_putfloat(Fdat, "d2", acpar->dr);
		sf_putfloat(Fdat, "o2", acpar->r0);
		sf_putstring(Fdat, "label2", "Receiver");
		sf_putstring(Fdat, "unit2", "km");
		sf_putint(Fdat, "n3", acpar->ns);
		sf_putfloat(Fdat, "d3", acpar->ds);
		sf_putfloat(Fdat, "o3", acpar->s0);
		sf_putstring(Fdat, "label3", "Shot");
		sf_putstring(Fdat, "unit3", "km");

		if(media==1) forward_modeling_a(Fdat, &mpipar, soupar, acpar, array, verb);
		else forward_modeling(Fdat, &mpipar, soupar, acpar, array, verb);

		sf_fileclose(Fdat);
	}
	else if(function == 2){ // FWI

		fwipar=(sf_fwi)sf_alloc(1, sizeof(*fwipar));
		if(!sf_getbool("onlygrad", &fwipar->onlygrad)) fwipar->onlygrad=false; 
		/* only calculate gradident or not */
		if(!sf_getint("grad_type",&fwipar->grad_type)) fwipar->grad_type=1;
		/* if 1, velocity; if 2, Q; if 3, velocity and Q */
		if(!sf_getint("prec",&fwipar->prec)) fwipar->prec=0;
		/* if 1, precondition gradient by dividing gradient locally by the energy of incident wavefield */
		if(!sf_getint("misfit_type", &fwipar->misfit_type)) fwipar->misfit_type=1;
		/* if 1, conventional misfit; if 2, wiener filter misfit; if 3, adaptive matching filtering misfit; if 4, smoothed misfit function */
		if(!sf_getint("match", &fwipar->match)) fwipar->match=1;
		/* when misfit_type=3, if match=0, use time-shift technique; if match=1, stretching */
		if(!sf_getbool("match_opt", &fwipar->match_opt)) fwipar->match_opt=false;
		/* apply optimization to get adjoint source */
		if(!sf_getint("fniter", &fwipar->fniter)) fwipar->fniter=100;
		/* iteration number of obtaining matching filter */
		if(!sf_getint("frectx", &fwipar->frectx)) fwipar->frectx=30;
		/* smoothing regularization radius in x */
		if(!sf_getint("frectz", &fwipar->frectz)) fwipar->frectz=50;
		/* smoothing regularization radius in z */
		if(!sf_getint("nw", &fwipar->nw)) fwipar->nw=21;
		/* number of time-shifts or stretching */
		if(!sf_getint("dw", &fwipar->dw)) fwipar->dw=10;
		/* interval of time-shifts or stretching */
		if(!sf_getint("drectx", &fwipar->drectx)) fwipar->drectx=1;
		/* smoothing radius in x for data residual */
		if(!sf_getint("drectz", &fwipar->drectz)) fwipar->drectz=50;
		/* smoothing radius in z for data residual */

		if(!sf_getfloat("wt1", &fwipar->wt1)) fwipar->wt1=acpar->t0;
		if(!sf_getfloat("wt2", &fwipar->wt2)) fwipar->wt2=acpar->t0+(acpar->nt-1)*acpar->dt;
		if(!sf_getfloat("woff1", &fwipar->woff1)) fwipar->woff1=acpar->r0;
		if(!sf_getfloat("woff2", &fwipar->woff2)) fwipar->woff2=acpar->r0+(acpar->nr-1)*acpar->dr;
		if(!sf_getfloat("gain", &fwipar->gain)) fwipar->gain=1;
		if(!sf_getint("waterz", &fwipar->waterz)) fwipar->waterz=51; /* water layer depth */
		if(!sf_getint("grectx", &fwipar->rectx)) fwipar->rectx=3; /* gradient smoothing radius in x */
		if(!sf_getint("grectz", &fwipar->rectz)) fwipar->rectz=3; /* gradient smoothing radius in z */

		Fdat=sf_input("Fdat"); /* input data */
		if(!fwipar->onlygrad){
			Finv=sf_output("output"); /* FWI result */
			Ferr=sf_output("Ferr"); /* data misfit convergence curve */
		}
		Fgrad=sf_output("Fgrad"); /* FWI gradient at first iteration */

		sf_putint(Fgrad, "n1", acpar->nz);
		sf_putfloat(Fgrad, "d1", acpar->dz);
		sf_putfloat(Fgrad, "o1", acpar->z0);
		sf_putstring(Fgrad, "label1", "Depth");
		sf_putstring(Fgrad, "unit1", "km");
		sf_putint(Fgrad, "n2", acpar->nx);
		sf_putfloat(Fgrad, "d2", acpar->dx);
		sf_putfloat(Fgrad, "o2", acpar->x0);
		sf_putstring(Fgrad, "label2", "Distance");
		sf_putstring(Fgrad, "unit2", "km");
		if(fwipar->grad_type==3) sf_putint(Fgrad, "n3", 2);

		if(!fwipar->onlygrad){
			optpar=(sf_optim)sf_alloc(1, sizeof(*optpar));
			if(!sf_getint("niter", &optpar->niter)) sf_error("iteration number required"); /* iteration number */
			if(!sf_getfloat("conv_error", &optpar->conv_error)) sf_error("convergence error required"); /* final convergence error */
			if(!sf_getint("opt_type", &optpar->opt_type)) optpar->opt_type=1;
			/* if 1, l-BFGS; if 2, steepest descent; if 3, multistep gradient */
			if(!sf_getint("npair", &optpar->npair)) optpar->npair=20; /* number of l-BFGS pairs */
			if(!sf_getint("nls", &optpar->nls)) optpar->nls=20; /* line search number */
			if(!sf_getfloat("factor", &optpar->factor)) optpar->factor=10; /* step length increase factor */
			if(!sf_getint("repeat", &optpar->repeat)) optpar->repeat=5; /* after how many iterations the step length goes back to 1 */
			if(!sf_getfloat("v1", &optpar->v1)) optpar->v1=0.;
			if(!sf_getfloat("v2", &optpar->v2)) optpar->v2=10.;
			if(!sf_getfloat("tau1", &optpar->tau1)) optpar->tau1=0.;
			if(!sf_getfloat("tau2", &optpar->tau2)) optpar->tau2=0.2;
			optpar->c1=1e-4;
			optpar->c2=0.9;
			if(!sf_getint("tangent", &optpar->tangent)) optpar->tangent=0; /* scheme 2 in smoothing kernel continuation */
			if(!sf_getint("sigma1", &optpar->sigma1)) optpar->sigma1=80; /* the change of 1st smoothing parameter */
			if(!sf_getint("sigma2", &optpar->sigma2)) optpar->sigma2=30; /* the change of 1st smoothing parameter */
			optpar->err=sf_floatalloc(optpar->niter+1);
		}
		/* dimension set up */
		if(Finv != NULL){
			sf_putint(Finv, "n1", acpar->nz);
			sf_putfloat(Finv, "d1", acpar->dz);
			sf_putfloat(Finv, "o1", acpar->z0);
			sf_putstring(Finv, "label1", "Depth");
			sf_putstring(Finv, "unit1", "km");
			sf_putint(Finv, "n2", acpar->nx);
			sf_putfloat(Finv, "d2", acpar->dx);
			sf_putfloat(Finv, "o2", acpar->x0);
			sf_putstring(Finv, "label2", "Distance");
			sf_putstring(Finv, "unit2", "km");
			if(fwipar->grad_type==3) sf_putint(Finv, "n3", 2);
			
			sf_putint(Ferr, "n1", optpar->niter+1);
			sf_putfloat(Ferr, "d1", 1);
			sf_putfloat(Ferr, "o1", 0);
			sf_putstring(Ferr, "label1", "Iterations");
			sf_putstring(Ferr, "unit1", "");
			sf_putint(Ferr, "n2", 1);
			sf_putfloat(Ferr, "d2", 1);
			sf_putfloat(Ferr, "o2", 0);
		}

		fwi(Fdat, Finv, Ferr, Fgrad, &mpipar, soupar, acpar, array, fwipar, optpar, verb, media);

		if(!fwipar->onlygrad){
			sf_fileclose(Finv);
			sf_fileclose(Ferr);
		}
		sf_fileclose(Fgrad);
	}
	else if(function == 3){ // RTM
		Fdat=sf_input("Fdat"); /* input data */
		Fimg=sf_output("output"); /* rtm image */

		/* dimension set up */
		sf_putint(Fimg, "n1", acpar->nz);
		sf_putfloat(Fimg, "d1", acpar->dz);
		sf_putfloat(Fimg, "o1", acpar->z0);
		sf_putstring(Fimg, "label1", "Depth");
		sf_putstring(Fimg, "unit1", "km");
		sf_putint(Fimg, "n2", acpar->nx);
		sf_putfloat(Fimg, "d2", acpar->dx);
		sf_putfloat(Fimg, "o2", acpar->x0);
		sf_putstring(Fimg, "label2", "Distance");
		sf_putstring(Fimg, "unit2", "km");

		if(media==1) rtm_a(Fdat, Fimg, &mpipar, soupar, acpar, array, verb);
		else rtm(Fdat, Fimg, &mpipar, soupar, acpar, array, verb);

		sf_fileclose(Fimg);
	}

	MPI_Finalize();
	exit(0);
}
