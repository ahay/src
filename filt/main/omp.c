/* OpenMP wrapper for embarassingly parallel jobs. */
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

#include <omp.h>

#include <rsf.h>

#define CMDLEN 4096

int main(int argc, char* argv[])
{
    int rank, nodes, node,ndim,n[SF_MAX_DIM],last,extra,chunk,i,j,len,nc;
    off_t size, size2, left,nbuf;
    char command[CMDLEN], *iname, *oname, key[5];
    char **inames, **onames, **cmdline, buffer[BUFSIZ];
    FILE *ifile, *ofile;
    sf_file inp, out, in;

#pragma omp parallel
    nodes = omp_get_num_threads();
    if (1 >= nodes) nodes = omp_get_num_procs(); 
#pragma omp end parallel

    sf_warning("Running %d threads",nodes);
   
    /* master node */
    sf_init(argc,argv);

    inp = sf_input("in");
    out = sf_output("out");

    ndim = sf_filedims (inp,n);
    size = sf_esize(inp);
    for (i=0; i < ndim-1; i++) {
	size *= n[i];
    }

    last = n[ndim-1];
    chunk = last/nodes;
    extra = last-chunk*nodes;
    snprintf(key,5,"n%d",ndim);

    j=0;
    for (i=1; i < argc; i++) {
	len = strlen(argv[i]);
	if (j+len > CMDLEN-2) sf_error("command line is too long");
	strncpy(command+j,argv[i],len);
	command[j+len]=' ';
	j += len+1;
    }
    command[j]='\0';

    inames = (char**) sf_alloc(nodes,sizeof(char*));
    onames = (char**) sf_alloc(nodes,sizeof(char*));
    cmdline = (char**) sf_alloc(nodes,sizeof(char*));

    for (node=0; node < nodes; node++) {
	if (node > 0) {
	    fclose(ifile);
	    fclose(ofile);
	}

	ifile = sf_tempfile(&iname,"w+b");
	ofile = sf_tempfile(&oname,"w+b");

	inames[node]=iname;
	onames[node]=oname;

	in = sf_output(iname);
	nc = (node==nodes-1)? chunk+extra:chunk;
	sf_putint(in,key,nc);
	sf_fileflush(in,inp);
	sf_setform(in,SF_NATIVE);

	for (i=0; i < nc; i++) {
	    for (nbuf=BUFSIZ,left=size; left > 0; left -= nbuf) {
		if (nbuf > left) nbuf=left;
		
		sf_charread(buffer,nbuf,inp);
		sf_charwrite(buffer,nbuf,in);
	    }
	}
	    
	sf_fileclose(in);

	cmdline[node] = sf_charalloc(CMDLEN);
	snprintf(cmdline[node],CMDLEN,"%s < %s > %s",command,iname,oname);
    }
    
#pragma omp parallel private(rank,sys)
    omp_set_num_threads(nodes);
    rank = omp_get_thread_num();
    sf_system(cmdline[rank]);
#pragma omp end parallel

    ifile = sf_tempfile(&iname,"w+b");
    ofile = sf_tempfile(&oname,"w+b");

    in = sf_output(iname);
    sf_cp(inp,in);
    sf_fileclose(inp);

    snprintf(command,CMDLEN,"%s dryrun=y < %s > %s",command,iname,oname);
    sf_system(command);
    sf_rm(iname,true,false,false);
    
    inp = sf_input(oname);
    ndim = sf_filedims (inp,n);
    if (last != n[ndim-1]) sf_error("Wrong dimensionality %d != %d",last,n[ndim-1]);

    size = sf_esize(inp);
    for (i=0; i < ndim-1; i++) {
	size *= n[i];
    }

    sf_setformat(out,sf_histstring(inp,"data_format"));
    sf_fileflush(out,inp);
    sf_setform(out,SF_NATIVE);
    sf_fileclose(inp);
    sf_rm(oname,true,false,false);

    for (node=0; node < nodes; node++) {
	iname = inames[node];
	oname = onames[node];
	    
	in = sf_input(oname);
	sf_setform(in,SF_NATIVE);

	nc = (node==nodes-1)? chunk+extra:chunk;
	for (i=0; i < nc; i++) {
	    for (nbuf=BUFSIZ,left=size; left > 0; left -= nbuf) {
		if (nbuf > left) nbuf=left;
		
		sf_charread(buffer,nbuf,in);
		sf_charwrite(buffer,nbuf,out);
	    }
	}

	sf_rm(iname,true,false,false);
	sf_rm(oname,true,false,false);
    }

    exit(0);
}
