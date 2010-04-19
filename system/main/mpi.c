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

#include <mpi.h>

#include <rsf.h>

#define CMDLEN 4096

static void sizes(sf_file file, int axis, int ndim, const off_t *n, off_t *size1, off_t *size2)
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
    int rank, nodes, skip, ndim,split,job,jobs,bigjobs,chunk,w,i,j,len, axis, axis2, splitargc;
    off_t size1, i2, size2, left,nbuf, n[SF_MAX_DIM], *splitsize1=NULL, *splitsize2=NULL;
    char cmdline[CMDLEN], command[CMDLEN], splitcommand[CMDLEN], *iname=NULL, *oname=NULL, key[5], *splitname;
    char **inames=NULL, **onames=NULL, buffer[BUFSIZ], **splitargv, *eq, *arg, *filekey;
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
	if (axis > ndim) axis=ndim;
	snprintf(key,5,"n%d",axis);

	axis--;
	sizes(inp,axis,ndim,n,&size1,&size2);
	split = n[axis];

	jobs = nodes-1;
	if (nodes < split) {
	    w = 1+(float) split/jobs;
	} else {
	    w = 1;
	    jobs = split;
	}
	bigjobs = split - jobs*(w-1);
	
	splitargv = (char**) sf_alloc(argc,sizeof(char*));
	splitargc = 0;

	j=0;
	for (i=1; i < argc; i++) {
	    if (strncmp(argv[i],"input=",6) &&
		strncmp(argv[i],"output=",7) &&
		strncmp(argv[i],"split=",6) &&
		strncmp(argv[i],"join=",5)) {
		if ('_'==argv[i][0]) {
		    splitargv[splitargc] = argv[i];
		    splitargc++;
		} else {
		    len = strlen(argv[i]);
		    if (j+len > CMDLEN-2) {
			sf_warning("command line is too long");
			MPI_Finalize();
		    }
		    strncpy(command+j,argv[i],len);
		    command[j+len]=' ';
		    j += len+1;
		}
	    }
	}
	command[j]='\0';

	if (0 < splitargc) { /* files to split other than input */
	    splitsize1 = sf_largeintalloc(splitargc);
	    splitsize2 = sf_largeintalloc(splitargc);
	    splitfile = (sf_file*) sf_alloc(splitargc,sizeof(sf_file));

	    for (i=0; i < splitargc; i++) {
		arg = splitargv[i];
		eq  = strchr(arg,'=');
		if (NULL == eq) {
		    sf_warning("Wrong parameter \"_%s\"",arg);
		    MPI_Finalize();
		}
		len = eq-arg;
		filekey = sf_charalloc(len+1);
		strncpy(filekey,arg,len);
		filekey[len]='\0';
	 
		splitfile[i] = sf_input(filekey);

		ndim = sf_largefiledims (splitfile[i],n);
		sizes(splitfile[i],axis,ndim,n,splitsize1+i,splitsize2+i);
	    }
	}

	inames = (char**) sf_alloc(jobs,sizeof(char*));
	onames = (char**) sf_alloc(jobs,sizeof(char*));

	for (job=0; job < jobs; job++) {
	    if (job < bigjobs) {
		chunk = w;
		skip = job*w;
	    } else {
		chunk = w-1;
		skip = bigjobs*w+(job-bigjobs)*chunk;
	    }

	    ifile = sf_tempfile(&iname,"w+b");
	    ofile = sf_tempfile(&oname,"w+b");

	    inames[job]=iname;
	    onames[job]=oname;

	    in = sf_output(iname);
	    fclose(ifile);

	    sf_putint(in,key,chunk);
	    sf_fileflush(in,inp);
	    sf_setform(in,SF_NATIVE);

	    for (i2=0; i2 < size2; i2++) {
		sf_seek(inp,size1*(i2*split+skip),SEEK_SET);

		nbuf=BUFSIZ;
		for (left=chunk*size1; left > 0; left -= nbuf) {
		    if (nbuf > left) nbuf=left;
		    
		    sf_charread(buffer,nbuf,inp);
		    sf_charwrite(buffer,nbuf,in);
		}
	    }

	    sf_fileclose(in);

	    splitcommand[0]='\0';
	    for (i=0; i < splitargc; i++) {
		arg = splitargv[i];
		eq  = strchr(arg,'=');
		len = eq-arg;
		filekey = sf_charalloc(len+1);
		strncpy(filekey,arg,len);
		filekey[len]='\0';
		
		ifile = sf_tempfile(&splitname,"w+b");
		in = sf_output(splitname);
		fclose(ifile);

		sf_putint(in,key,chunk);
		sf_fileflush(in,splitfile[i]);
		sf_setform(in,SF_NATIVE);

		for (i2=0; i2 < size2; i2++) {
		    sf_seek(splitfile[i],size1*(i2*split+skip),SEEK_SET);

		    nbuf=BUFSIZ;
		    for (left=chunk*size1; left > 0; left -= nbuf) {
			if (nbuf > left) nbuf=left;
			
			sf_charread(buffer,nbuf,splitfile[i]);
			sf_charwrite(buffer,nbuf,in);
		    }
		}

		sf_fileclose(in);

		snprintf(cmdline,CMDLEN,"%s %s=%s",splitcommand,filekey+1,splitname);
		strncpy(splitcommand,cmdline,CMDLEN);
	    }	

	    snprintf(cmdline,CMDLEN,"%s %s < %s > %s",command,splitcommand,iname,oname);
	    MPI_Send(cmdline,strlen(cmdline)+1,MPI_CHAR,job+1,0,MPI_COMM_WORLD);
	    fclose(ofile);
	}

	iname = sf_getstring("input");
	ofile = sf_tempfile(&oname,"w+b");

	snprintf(cmdline,CMDLEN,"%s --dryrun=y < %s > %s",command,iname,oname);
	sf_system(cmdline);

	inp = sf_input(oname);
	ndim = sf_largefiledims (inp,n);

	if (!sf_getint("join",&axis2)) axis2=axis+1;
	if (axis2 > ndim) axis2=ndim;
	snprintf(key,5,"n%d",axis2);

	axis2--;
	sizes(inp,axis,ndim,n,&size1,&size2);

	sf_setformat(out,sf_histstring(inp,"data_format"));
	sf_fileflush(out,inp);

	sf_setform(out,SF_NATIVE);
	sf_rm(oname,true,false,false);

	for (job=0; job < jobs; job++) {
	    MPI_Recv(&rank,1, MPI_INT, job+1, 1, MPI_COMM_WORLD,&stat);

	    iname = inames[job];
	    oname = onames[job];

	    in = sf_input(oname);
	    sf_setform(in,SF_NATIVE);

	    if (!sf_histint(in,key,&chunk)) {
		sf_warning("Problem on node %d with file %s",job+1,oname);
		MPI_Finalize();
	    }
	    
	    for (i2=0; i2 < size2; i2++) {
		nbuf=BUFSIZ;
		for (left=chunk*size1; left > 0; left -= nbuf) {
		    if (nbuf > left) nbuf=left;
		    
		    sf_charread(buffer,nbuf,in);
		    sf_charwrite(buffer,nbuf,out);
		}
	    }

	    sf_fileclose(in);

	    sf_rm(iname,true,false,false);
	    sf_rm(oname,true,false,false);
	}

	sf_fileclose(inp);
    } else { /* slave nodes */
	MPI_Recv(cmdline, CMDLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&stat);
	sf_warning(cmdline);
        sf_system(cmdline);
	MPI_Send(&rank,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
}
