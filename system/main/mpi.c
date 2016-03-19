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

int main(int argc, char* argv[])
{
    int rank, nodes, ndim, job, axis, axis2, jobs;
    off_t n[SF_MAX_DIM];
    char **commands, cmdline[SF_CMDLEN], *iname;
    sf_file inp=NULL, out=NULL;
    MPI_Status stat;

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    if (nodes < 2) {
	fprintf(stderr,"Need at least two nodes!\n");
	MPI_Finalize();
	exit(1);
    }

    if (!rank) { /* master node */
	sf_init(argc,argv);

	inp = sf_input("--input");
	out = sf_output("--output");

	ndim = sf_largefiledims (inp,n);

	if (!sf_getint("split",&axis)) axis=ndim;
	/* axis to split */
	
	commands = sf_split(inp,axis,nodes,&jobs,ndim,n,argc,argv);  

	for (job=0; job < jobs; job++) {
	    strncpy(cmdline,commands[job],SF_CMDLEN);
	    MPI_Send(cmdline, SF_CMDLEN, MPI_CHAR, job+1, 0, MPI_COMM_WORLD);
	}

	iname = sf_getstring("--input");

	if (!sf_getint("join",&axis2)) axis2=axis;
	/* axis to join (0 means add) */

	sf_out(out,jobs,axis2,iname);
	
	for (job=0; job < jobs; job++) {
	    MPI_Recv(&rank,1, MPI_INT, job+1, 1, MPI_COMM_WORLD,&stat);
	}
	
	if (axis2 > 0) {
	    sf_join(out,axis2,jobs);
	} else {
	    sf_add(out,jobs);
	}

	sf_fileclose(inp);
    } else { /* slave nodes */
	MPI_Recv(cmdline, SF_CMDLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&stat);
	fprintf(stderr,"CPU %d: %s\n",rank,cmdline); 
        sf_system(cmdline);
	MPI_Send(&rank,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
}
