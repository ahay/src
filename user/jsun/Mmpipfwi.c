/* Visco-acoustic Forward Modeling, FWI, and RTM based on SLS model */
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Qfwi_modeling.h"
#include "Qfwi_fwi.h"
#include "Qfwi_rtm.h"

int main(int argc, char* argv[])
{
	int function, media, ntmp;
	bool verb;

	sf_mpi mpipar;
	sf_sou soupar;
	sf_acqui acpar;
	sf_vec array;
	sf_fwi fwipar = NULL;
	sf_optim optpar=NULL;
	sf_pas paspar=NULL;

	MPI_Comm comm=MPI_COMM_WORLD;

	sf_file Fv, Fq, Fw, Fdat, Fimg, Finv=NULL, Fgrad=NULL, Fsrc, Fmwt=NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &mpipar.cpuid);
	MPI_Comm_size(comm, &mpipar.numprocs);

	sf_init(argc, argv);

#ifdef _OPENMP
        omp_init();
#endif

	Fv=sf_input("Fvel"); /* velocity model */
	Fq=sf_input("Fq"); /* quality factor */
	Fw=sf_input("Fwavelet"); /* wavelet */

	soupar=(sf_sou)sf_alloc(1, sizeof(*soupar));
	acpar=(sf_acqui)sf_alloc(1, sizeof(*acpar));
	array=(sf_vec)sf_alloc(1, sizeof(*array));

	/* parameters I/O */
	if(!sf_getint("media", &media)) media=1;
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
	if(!sf_getfloat("coef", &acpar->coef)) acpar->coef=0.003; /* absorbing boundary coefficient */

	if(!sf_getint("acqui_type", &acpar->acqui_type)) acpar->acqui_type=1;
	/* if 1, fixed acquisition; if 2, marine acquisition; if 3, symmetric acquisition */
	if(!sf_getint("ns", &acpar->ns)) sf_error("shot number required"); /* shot number */
	if(!sf_getfloat("ds", &acpar->ds)) sf_error("shot interval required"); /* shot interval */
	if(!sf_getfloat("s0", &acpar->s0)) sf_error("shot origin required"); /* shot origin */
	if(!sf_getint("sz", &acpar->sz)) acpar->sz=5; /* source depth */
	if(!sf_getint("nr", &acpar->nr)) acpar->nr=acpar->nx; /* number of receiver */
	if(!sf_getfloat("dr", &acpar->dr)) acpar->dr=acpar->dx; /* receiver interval */
	if(!sf_getfloat("r0", &acpar->r0)) acpar->r0=acpar->x0; /* receiver origin */
	if(!sf_getint("rz", &acpar->rz)) acpar->rz=1; /* receiver depth */

	if(!sf_getfloat("f0", &acpar->f0)) sf_error("reference frequency required"); /* reference frequency */
	if(!sf_getint("interval", &acpar->interval)) acpar->interval=1; /* wavefield storing interval */

	if(!sf_getfloat("fhi", &soupar->fhi)) soupar->fhi=0.5/acpar->dt; 
	if(!sf_getfloat("flo", &soupar->flo)) soupar->flo=0.; 
	soupar->rectx=2; 
	soupar->rectz=2; 

	/* get prepared */
	preparation(Fv, Fq, Fw, acpar, soupar, array);

        switch (function) {

            case 1: /* Modeling */

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

                break;

            case 2: /* FWI */

		fwipar=(sf_fwi)sf_alloc(1, sizeof(*fwipar));
		if(!sf_getbool("onlygrad", &fwipar->onlygrad)) fwipar->onlygrad=false; /* only want gradident */
		fwipar->grad_type=1;
		fwipar->misfit_type=1;
		fwipar->opt_type=1;
                if(!sf_getfloat("wt1", &fwipar->wt1)) fwipar->wt1=acpar->t0;
                if(!sf_getfloat("wt2", &fwipar->wt2)) fwipar->wt2=acpar->t0+(acpar->nt-1)*acpar->dt;
                if(!sf_getfloat("woff1", &fwipar->woff1)) fwipar->woff1=acpar->r0;
                if(!sf_getfloat("woff2", &fwipar->woff2)) fwipar->woff2=acpar->r0+(acpar->nr-1)*acpar->dr;
                if(!sf_getbool("oreo", &fwipar->oreo)) fwipar->oreo=false; /* keep oreo or keep cream */
		if(!sf_getint("waterz", &fwipar->waterz)) fwipar->waterz=51; /* water layer depth */
		if(!sf_getint("grectx", &fwipar->rectx)) fwipar->rectx=3; /* gradient smoothing radius in x */
		if(!sf_getint("grectz", &fwipar->rectz)) fwipar->rectz=3; /* gradient smoothing radius in z */

		Fdat=sf_input("Fdat"); /* input data */
		if(!fwipar->onlygrad) Finv=sf_output("output"); /* FWI result */
		Fgrad=sf_output("Fgrad"); /* FWI gradient at first iteration */

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
		}
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
			optpar->npair=20; /* number of l-BFGS pairs */
			optpar->nls=20; /* line search number */
			if(!sf_getfloat("c1", &optpar->c1)) optpar->c1=1e-4;
			if(!sf_getfloat("c2", &optpar->c2)) optpar->c2=0.9;
			optpar->factor=10;
                        if(!sf_getfloat("v1", &optpar->v1)) optpar->v1=0.;
                        if(!sf_getfloat("v2", &optpar->v2)) optpar->v2=10.;
		}

		fwi(Fdat, Finv, Fgrad, &mpipar, soupar, acpar, array, fwipar, optpar, verb, media);

		if(!fwipar->onlygrad) sf_fileclose(Finv);
		sf_fileclose(Fgrad);

                break;

            case 3: /* RTM */

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
                
                break;

            case 4: /* Passive FWI */

                paspar = passive_init(acpar);

                if (paspar->inv) {

                    Fdat=sf_input("Fdat");
                    if(!sf_histint(Fdat, "n3", &acpar->ns)) acpar->ns=1;

                    if (!paspar->onlysrc) {
                        fwipar=(sf_fwi)sf_alloc(1, sizeof(*fwipar));
                        if(!sf_getbool("onlygrad", &fwipar->onlygrad)) fwipar->onlygrad=false; /* only want gradident */
                        fwipar->grad_type=1;
                        fwipar->misfit_type=1;
                        fwipar->opt_type=1;
                        if(!sf_getfloat("wt1", &fwipar->wt1)) fwipar->wt1=acpar->t0;
                        if(!sf_getfloat("wt2", &fwipar->wt2)) fwipar->wt2=acpar->t0+(acpar->nt-1)*acpar->dt;
                        if(!sf_getfloat("woff1", &fwipar->woff1)) fwipar->woff1=acpar->r0;
                        if(!sf_getfloat("woff2", &fwipar->woff2)) fwipar->woff2=acpar->r0+(acpar->nr-1)*acpar->dr;
                        if(!sf_getbool("oreo", &fwipar->oreo)) fwipar->oreo=false; /* keep oreo or keep cream */
                        if(!sf_getint("waterz", &fwipar->waterz)) fwipar->waterz=0; /* water layer depth */
                        if(!sf_getint("waterzb", &fwipar->waterzb)) fwipar->waterzb=0; /* water layer depth from bottom up */
                        if(!sf_getint("grectx", &fwipar->rectx)) fwipar->rectx=3; /* gradient smoothing radius in x */
                        if(!sf_getint("grectz", &fwipar->rectz)) fwipar->rectz=3; /* gradient smoothing radius in z */

                        if(!fwipar->onlygrad) Finv=sf_output("output"); /* FWI result */
                        Fgrad=sf_output("Fgrad"); /* FWI gradient at first iteration */

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
                            /*if(fwipar->grad_type==3) sf_putint(Finv, "n3", 2);*/
                        }
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
                        /*if(fwipar->grad_type==3) sf_putint(Fgrad, "n3", 2);*/

                        if(!fwipar->onlygrad){
                            optpar=(sf_optim)sf_alloc(1, sizeof(*optpar));
                            if(!sf_getint("niter", &optpar->niter)) sf_error("iteration number required"); /* iteration number */
                            if(!sf_getint("repeat", &optpar->repeat)) optpar->repeat=1; /* repeat resetting alpha */
                            if(!sf_getfloat("conv_error", &optpar->conv_error)) sf_error("convergence error required"); /* final convergence error */
                            optpar->npair=20; /* number of l-BFGS pairs */
                            optpar->nls=20; /* line search number */
                            if(!sf_getfloat("c1", &optpar->c1)) optpar->c1=1e-4;
                            if(!sf_getfloat("c2", &optpar->c2)) optpar->c2=0.9;
                            optpar->factor=10;
                            if(!sf_getfloat("v1", &optpar->v1)) optpar->v1=0.;
                            if(!sf_getfloat("v2", &optpar->v2)) optpar->v2=10.;
                        }
                    } /* if !onlysrc */

                    if (!paspar->onlyvel) {
                        Fsrc=sf_output("Fsrc");
                        sf_putint   (Fsrc, "n1", acpar->nz);
                        sf_putfloat (Fsrc, "o1", acpar->z0);
                        sf_putfloat (Fsrc, "d1", acpar->dz);
                        sf_putstring(Fsrc, "label1", "Depth");
                        sf_putstring(Fsrc, "unit1" , "km");
                        sf_putint   (Fsrc, "n2", acpar->nx);
                        sf_putfloat (Fsrc, "o2", acpar->x0);
                        sf_putfloat (Fsrc, "d2", acpar->dx);
                        sf_putstring(Fsrc, "label2", "Distance");
                        sf_putstring(Fsrc, "unit2" , "km");
                        sf_putint   (Fsrc, "n3", acpar->nt);
                        sf_putfloat (Fsrc, "o3", acpar->t0);
                        sf_putfloat (Fsrc, "d3", acpar->dt);
                        sf_putstring(Fsrc, "label3", "Time");
                        sf_putstring(Fsrc, "unit3" , "s");
                        sf_putint   (Fsrc, "n4", acpar->ns);
                        sf_putfloat (Fsrc, "d4", 1.0f);
                        sf_putfloat (Fsrc, "o4", 0.0f);
                        sf_putstring(Fsrc, "label4", "Stage");

                        Fmwt=sf_output("Fmwt"); /* output data */
                        sf_putint   (Fmwt, "n1", acpar->nz);
                        sf_putfloat (Fmwt, "o1", acpar->z0);
                        sf_putfloat (Fmwt, "d1", acpar->dz);
                        sf_putstring(Fmwt, "label1", "Depth");
                        sf_putstring(Fmwt, "unit1" , "km");
                        sf_putint   (Fmwt, "n2", acpar->nx);
                        sf_putfloat (Fmwt, "o2", acpar->x0);
                        sf_putfloat (Fmwt, "d2", acpar->dx);
                        sf_putstring(Fmwt, "label2", "Distance");
                        sf_putstring(Fmwt, "unit2" , "km");
                        sf_putint   (Fmwt, "n3", acpar->nt);
                        sf_putfloat (Fmwt, "o3", acpar->t0);
                        sf_putfloat (Fmwt, "d3", acpar->dt);
                        sf_putstring(Fmwt, "label3", "Time");
                        sf_putstring(Fmwt, "unit3" , "s");
                    } else {
                        Fsrc=sf_input("Fsrc");
                        if(!sf_histint(Fsrc, "n4", &ntmp)) ntmp=1;
                        if (ntmp!=acpar->ns) sf_error("Shot dimension mismatch!");
                    }

                } else { /* modeling */
                    Fsrc=sf_input("Fsrc");
                    if(!sf_histint(Fsrc, "n4", &acpar->ns)) acpar->ns=1;

                    Fdat=sf_output("output"); /* output data */
                    sf_putint   (Fdat, "n1", acpar->nt);
                    sf_putfloat (Fdat, "o1", acpar->t0);
                    sf_putfloat (Fdat, "d1", acpar->dt);
                    sf_putstring(Fdat, "label1", "Time");
                    sf_putstring(Fdat, "unit1" , "s");
                    sf_putint   (Fdat, "n2", acpar->nx);
                    sf_putfloat (Fdat, "o2", acpar->x0);
                    sf_putfloat (Fdat, "d2", acpar->dx);
                    sf_putstring(Fdat, "label2", "Distance");
                    sf_putstring(Fdat, "unit2" , "km");
                    sf_putint   (Fdat, "n3", acpar->ns);
                    sf_putfloat (Fdat, "d3", 1.0f);
                    sf_putfloat (Fdat, "o3", 0.0f);
                    sf_putstring(Fdat, "label3", "Stage");
                }

                if (paspar->inv) {
                    if (paspar->onlysrc) { /* only inverting for source */
                        lstri(Fdat, Fmwt, Fsrc, &mpipar, acpar, array, paspar, verb);
                        sf_fileclose(Fsrc);
                        sf_fileclose(Fmwt);
                    } else { /* inverting for velocity ( and source ) */
                        pfwi(Fdat, Finv, Fgrad, Fmwt, Fsrc, &mpipar, soupar, acpar, array, fwipar, optpar, paspar, verb);
                        if(!fwipar->onlygrad) sf_fileclose(Finv);
                        sf_fileclose(Fgrad);
                        if (!paspar->onlyvel) {
                            sf_fileclose(Fsrc); 
                            sf_fileclose(Fmwt);
                        }
                    }
                } else {
                    lstri(Fdat, Fmwt, Fsrc, &mpipar, acpar, array, paspar, verb);
                    sf_fileclose(Fdat);
                }

                break;

            default:
                sf_warning("Please specify a valid function");

	} /* switch */

	MPI_Finalize();
	exit(0);
}
