/* 2-D dip estimation by plane wave destruction with MPI parallelization. */
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
#include <rsfpwd.h>

#include <mpi.h>

int main (int argc, char *argv[])
{
    int n123, niter, order, nj1,nj2, i, j, liter, dim;
    int n[SF_MAX_DIM], rect[3], nr, ir; 
    float p0, *u, *p, pmin, pmax, eps;
    float **allu, **allp, d1, d2, d3, o1, o2, o3, *sendbuf, *recvbuf;
    bool verb, **mm;
    sf_file in, out, mask, dip0;

    int cpuid, numprocs, nrpad, iturn;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &cpuid);
    MPI_Comm_size(comm, &numprocs);

    sf_init(argc,argv);
    in = sf_input ("input");
    out = sf_output ("output");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");
    if (cpuid==0) sf_warning("numprocs=%d", numprocs);

    dim = sf_filedims(in,n);
    if (dim < 2) n[1]=1;
    n123 = n[0]*n[1];
    nr = 1;
    for (j=2; j < dim; j++) {
        nr *= n[j];
    }

    n[2]= 1;
    rect[2]=1;
    nj2=1;

    if(nr%numprocs==0)
        nrpad=nr;
    else
        nrpad=(nr/numprocs+1)*numprocs;

    /* set up output dimension */
    if(cpuid==0){
        if (!sf_histfloat(in, "d1", &d1)) sf_error("No d1= in input");
        if (!sf_histfloat(in, "o1", &o1)) sf_error("No o1= in input");
        if (!sf_histfloat(in, "d2", &d2)) sf_error("No d2= in input");
        if (!sf_histfloat(in, "o2", &o2)) sf_error("No o2= in input");
        if (!sf_histfloat(in, "d3", &d3)) sf_error("No d3= in input");
        if (!sf_histfloat(in, "o3", &o3)) sf_error("No o3= in input");

        sf_putint(out, "n1", n[0]);
        sf_putfloat(out, "d1", d1);
        sf_putfloat(out, "o1", o1);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "m");
        sf_putint(out, "n2", n[1]);
        sf_putfloat(out, "d2", d2);
        sf_putfloat(out, "o2", o2);
        sf_putstring(out, "label2", "Offset");
        sf_putstring(out, "unit2", "m");
        sf_putint(out, "n3", nr);
        sf_putfloat(out, "d3", d3);
        sf_putfloat(out, "o3", o3);
        sf_putstring(out, "label3", "CIGs");
        sf_putstring(out, "unit3", "m");
    }

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* dip smoothness on 2nd axis */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial dip */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("pmin",&pmin)) pmin = -FLT_MAX;
    /* minimum dip */
    if (!sf_getfloat("pmax",&pmax)) pmax = +FLT_MAX;
    /* maximum dip */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    /* initialize dip estimation */
    dip3_init(n[0], n[1], n[2], rect, liter, eps, false);

    /* initial dip file */
    if(NULL != sf_getstring("dip0")){
        dip0=sf_input("dip0");
    }else{
        dip0=NULL;
    }
    if(cpuid==0){
        allu=sf_floatalloc2(n123, nrpad);
        sf_floatread(allu[0], n123*nr, in);
        for(ir=nr; ir<nrpad; ir++)
            for(i=0; i<n123; i++)
                allu[ir][i]=0.;
        allp=sf_floatalloc2(n123, nrpad);
        if(NULL != dip0){
            sf_floatread(allp[0], n123*nr, dip0);
            for(ir=nr; ir<nrpad; ir++)
                for(i=0; i<n123; i++)
                    allp[ir][i]=0.;
        }
    }
    u = sf_floatalloc(n123);
    p = sf_floatalloc(n123);

    /* masking operator */
    if(NULL != sf_getstring("mask")) {
        mm = sf_boolalloc2(n123,2);
        mask = sf_input("mask");
        if(cpuid==0) sf_floatread(u, n123, mask);
        MPI_Bcast(u, n123, MPI_FLOAT, 0, comm);
        mask32 (false, order, nj1, nj2, n[0], n[1], n[2], u, mm);
    }else{
        mm = (bool**) sf_alloc(2,sizeof(bool*));
        mm[0] = mm[1] = NULL;
    }

    /* loop over third dimension */
    for(iturn=0; iturn*numprocs<nrpad; iturn++){
        ir=iturn*numprocs+cpuid;
        if (cpuid==0 && verb) sf_warning("slice %d of %d;", ir+1, nr);

        /* image data */
        if(cpuid==0){
            sendbuf=allu[iturn*numprocs];
            recvbuf=u;
        }else{
            sendbuf=NULL;
            recvbuf=u;
        }
        MPI_Scatter(sendbuf, n123, MPI_FLOAT, recvbuf, n123, MPI_FLOAT, 0, comm);

        /* initialize t-x dip */
        if(NULL != dip0) {
            if(cpuid==0){
                sendbuf=allp[iturn*numprocs];
                recvbuf=p;
            }else{
                sendbuf=NULL;
                recvbuf=p;
            }
            MPI_Scatter(sendbuf, n123, MPI_FLOAT, recvbuf, n123, MPI_FLOAT, 0, comm);
        }else{
            for(i=0; i < n123; i++) {
                p[i] = p0;
            }
        }

        /* estimate t-x dip */
        if(ir<nr) dip3(false, 1, niter, order, nj1, u, p, mm[0], pmin, pmax);

        if(cpuid==0){
            sendbuf=p;
            recvbuf=allp[iturn*numprocs];
        }else{
            sendbuf=p;
            recvbuf=NULL;
        }
        MPI_Gather(sendbuf, n123, MPI_FLOAT, recvbuf, n123, MPI_FLOAT, 0, comm);
    }

    if(cpuid==0) sf_floatwrite(allp[0], n123*nr, out);

    MPI_Finalize();
    exit (0);
}
