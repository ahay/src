/* 2-D painting by plane-wave construction with MPI parallelization. */
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
    bool verb;
    int n1,n2,n3, i1,i2,i3, i0, order;
    float o1, d1, eps, **u, **p, *trace, *time;
    float ***allu, ***allp, d2, d3, o2, o3, *sendbuf, *recvbuf;
    sf_file out, dip;

    int cpuid, numprocs, nrpad, iturn, n12;
    MPI_Comm comm=MPI_COMM_WORLD;

    sf_init(argc,argv);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &cpuid);
    MPI_Comm_size(comm, &numprocs);

    dip = sf_input("--input");
    out = sf_output("--output");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(dip,2);
    n12=n1*n2;

    if(n3%numprocs==0)
        nrpad=n3;
    else
        nrpad=(n3/numprocs+1)*numprocs;

    /* reference trace */
    time = sf_floatalloc(n1);
    if (!sf_histfloat(dip,"o1",&o1)) o1=0.; 
    if (!sf_histfloat(dip,"d1",&d1)) d1=1.; 
    for (i1=0; i1 < n1; i1++) {
        time[i1] = o1+i1*d1;
    }

    /* set up output dimension */
    if(cpuid==0){
        if (!sf_histfloat(dip, "d2", &d2)) sf_error("No d2= in input");
        if (!sf_histfloat(dip, "o2", &o2)) sf_error("No o2= in input");
        if (!sf_histfloat(dip, "d3", &d3)) sf_error("No d3= in input");
        if (!sf_histfloat(dip, "o3", &o3)) sf_error("No o3= in input");

        sf_putint(out, "n1", n1);
        sf_putfloat(out, "d1", d1);
        sf_putfloat(out, "o1", o1);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1", "m");
        sf_putint(out, "n2", n2);
        sf_putfloat(out, "d2", d2);
        sf_putfloat(out, "o2", o2);
        sf_putstring(out, "label2", "Offset");
        sf_putstring(out, "unit2", "m");
        sf_putint(out, "n3", n3);
        sf_putfloat(out, "d3", d3);
        sf_putfloat(out, "o3", o3);
        sf_putstring(out, "label3", "CIGs");
        sf_putstring(out, "unit3", "m");
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("i0",&i0)) i0=0;
    /* reference trace */
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1, n2, eps*eps, order, 1, false);

    if(cpuid==0){
        allp = sf_floatalloc3(n1, n2, nrpad);
        sf_floatread(allp[0][0], n12*n3, dip);
        for(i3=n3; i3<nrpad; i3++)
            for(i2=0; i2<n2; i2++)
                for(i1=0; i1<n1; i1++)
                    allp[i3][i2][i1]=0.;
        allu = sf_floatalloc3(n1, n2, nrpad);
    }
    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n1);

    for(iturn=0; iturn*numprocs<nrpad; iturn++){
        i3=iturn*numprocs+cpuid;
        if(cpuid==0 && verb) sf_warning("slice %d of %d;", i3+1, n3);

        /* scatter data from master to other processors */
        if(cpuid==0){
            sendbuf=allp[iturn*numprocs][0];
            recvbuf=p[0];
        }else{
            sendbuf=NULL;
            recvbuf=p[0];
        }
        MPI_Scatter(sendbuf, n12, MPI_FLOAT, recvbuf, n12, MPI_FLOAT, 0, comm);

        if(i3<n3){ /* predictive painting */
            for (i1=0; i1 < n1; i1++) {
                trace[i1] = time[i1];
            }

            for (i1=0; i1 < n1; i1++) {
                u[i0][i1] = trace[i1];
            }
            for (i2=i0-1; i2 >= 0; i2--) {
                predict_step(false,false,trace,p[i2]);
                for (i1=0; i1 < n1; i1++) {
                    u[i2][i1] = trace[i1];
                }
            }
            for (i1=0; i1 < n1; i1++) {
                trace[i1] = u[i0][i1];
            }
            for (i2=i0+1; i2 < n2; i2++) {
                predict_step(false,true,trace,p[i2-1]);
                for (i1=0; i1 < n1; i1++) {
                    u[i2][i1] = trace[i1];
                }
            }
        }

        /* gather data from other processors to master */
        if(cpuid==0){
            sendbuf=u[0];
            recvbuf=allu[iturn*numprocs][0];
        }else{
            sendbuf=u[0];
            recvbuf=NULL;
        }
        MPI_Gather(sendbuf, n12, MPI_FLOAT, recvbuf, n12, MPI_FLOAT, 0, comm);
    }

    if(cpuid==0) sf_floatwrite(allu[0][0], n12*n3, out);

    MPI_Finalize();
    exit (0);
}

/* 	$Id: Mflat.c 1131 2005-04-20 18:19:10Z fomels $	 */
