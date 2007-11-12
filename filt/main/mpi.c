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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>

#include <mpi.h>

#include <rsf.h>

#define CMDLEN 4096

int main(int argc, char* argv[])
{
    int rank, nodes, node,ndim,n[SF_MAX_DIM],last,chunk,i,j,len,nc,sys;
    off_t size;
    char cmdline[CMDLEN], command[CMDLEN], *iname, *oname, *data, key[5];
    char **inames, **onames;
    sf_file inp, out, in;
    MPI_Status stat;

    MPI_Init(&argc,&argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    
    if (!rank) { /* master node */
	sf_init(argc,argv);

	inp = sf_input("input");
	out = sf_output("output");

	ndim = sf_filedims (inp,n);
	size = sf_esize(inp);
	for (i=0; i < ndim-1; i++) {
	    size *= n[i];
	}

	data = sf_charalloc(size);

	last = n[ndim-1];
	chunk = last/(nodes-1);
	last = last-chunk*(nodes-1);
	snprintf(key,5,"n%d",ndim);

	j=0;
	for (i=1; i < argc; i++) {
	    if (strncmp(argv[i],"input=",6) &&
		strncmp(argv[i],"output=",7)) {
		len = strlen(argv[i]);
		if (j+len > CMDLEN-2) sf_error("command line is too long");
		strncpy(command+j,argv[i],len);
		command[j+len]=' ';
		j += len+1;
	    }
	}
	command[j]='\0';

	inames = (char**) sf_alloc(nodes-1,sizeof(char*));
	onames = (char**) sf_alloc(nodes-1,sizeof(char*));

	for (node=1; node < nodes; node++) {
	    iname = sf_tempname();
	    oname = sf_tempname();

	    inames[node-1]=iname;
	    onames[node-1]=oname;

	    in = sf_output(iname);

	    nc = (node==nodes-1)? chunk+last:chunk;

	    sf_putint(in,key,nc);
	    sf_fileflush(in,inp);
	    sf_setform(in,SF_NATIVE);

	    for (i=0; i < nc; i++) {
		sf_charread(data,size,inp);
		sf_charwrite(data,size,in);
	    }
	    
	    sf_fileclose(in);

	    snprintf(cmdline,CMDLEN,"%s < %s > %s",command,iname,oname);
	    MPI_Send(cmdline,strlen(cmdline)+1,MPI_CHAR,node,0,MPI_COMM_WORLD);
	}

	sf_fileflush(out,inp);
	sf_setform(out,SF_NATIVE);

	for (node=1; node < nodes; node++) {
	    MPI_Recv(&sys,1, MPI_INT, node, 1, MPI_COMM_WORLD,&stat);

	    iname = inames[node-1];
	    oname = onames[node-1];
	    
	    in = sf_input(oname);
	    sf_setform(in,SF_NATIVE);

	    nc = (node==nodes-1)? chunk+last:chunk;
	    for (i=0; i < nc; i++) {
		sf_charread(data,size,in);
		sf_charwrite(data,size,out);
	    }

	    sf_rm(iname,true,false,false);
	    sf_rm(oname,true,false,false);
	}
    } else { /* slave nodes */
	MPI_Recv(cmdline, CMDLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&stat);
        sys = system(cmdline);
	if (sys ==-1) printf("\"%s\" failed on Node %d\n",cmdline,rank);
	MPI_Send(&sys,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
    
