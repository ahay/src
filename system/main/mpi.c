/* MPI wrapper for embarassingly parallel jobs. */
/*
  Copyright (C) 2007 University of Texas at Austin

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
#include <mpi.h>

#include <rsf.h>

#define CMDLEN 4096

static void sizes(sf_file file, int axis, int ndim, 
		  const off_t *n, off_t *size1, off_t *size2)
{
    int i;

    *size1 = sf_esize(file);
    *size2 = 1;
    for (i=0; i < ndim; i++) {
	if      (i < axis) *size1 *= n[i];
	else if (i > axis) *size2 *= n[i];
    }
}

int main(int argc, char* argv[])
{
    int rank, nodes, skip, ndim,split,job,jobs,bigjobs,chunk,w,i,j,len;
    int axis, axis2, splitargc;
    off_t size1, i2, size2, left,nbuf, n[SF_MAX_DIM];
    off_t *splitsize1=NULL, *splitsize2=NULL;
    float d, o, di, oi;
    char **commands, cmdline[CMDLEN], command[CMDLEN], splitcommand[CMDLEN];
    char *iname=NULL, *oname=NULL, nkey[5], okey[5], dkey[5], *splitname;
    char **inames=NULL, **onames=NULL, ***spnames=NULL;
    char buffer[BUFSIZ], **splitkey, *eq, *arg;
    FILE *ifile=NULL, *ofile=NULL;
    sf_file inp=NULL, out=NULL, in=NULL, *splitfile=NULL;
    MPI_Status stat;

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    if (nodes < 2) {
	fprintf(stderr,"Need at least two nodes!\n");
	MPI_Finalize();
    }

    if (!rank) { /* master node */
	sf_init(argc,argv);

	inp = sf_input("input");
	out = sf_output("output");

	ndim = sf_largefiledims (inp,n);

	if (!sf_getint("split",&axis)) axis=ndim;
	/* axis to split */
	
	commands = parallel_split(inp,axis,nodes,ndim,n,argc,argv);  

	for (job=1; job < nodes; job++) {
	    cmdline = commands[job-1];
	    MPI_Send(cmdline,strlen(cmdline)+1,MPI_CHAR,job,0,MPI_COMM_WORLD);
	}

	iname = sf_getstring("input");
	if (!sf_getint("join",&axis2)) axis2=axis;

	parallel_out(out,axis2,iname);
	
	for (job=1; job < nodes; job++) {
	    MPI_Recv(&rank,1, MPI_INT, job, 1, MPI_COMM_WORLD,&stat);
	    parallel_join(out,job-1);
	}

	sf_fileclose(inp);
    } else { /* slave nodes */
	MPI_Recv(cmdline, CMDLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&stat);
	fprintf(stderr,"CPU %d: %s\n",rank,cmdline); 
        sf_system(cmdline);
	MPI_Send(&rank,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
}
