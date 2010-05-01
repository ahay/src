/* Split parallel jobs. */

/*
  Copyright (C) 2010 University of Texas at Austin

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
#include <string.h>
#include <stdio.h>

#include "parallel.h"

#include "file.h"
/*^*/

#include "alloc.h"
#include "error.h"
#include "files.h"
#include "system.h"

#define SF_CMDLEN 4096
/*^*/

static char command[SF_CMDLEN], splitcommand[SF_CMDLEN], **inames=NULL, **onames=NULL;
static char ***spnames=NULL, buffer[BUFSIZ], nkey[5];
static off_t size1, size2;
static int splitargc;

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

char** sf_split(sf_file inp          /* input file */, 
		      int axis             /* split axis */,
		      int nodes            /* number of CPUs */,
		      int ndim, off_t *n   /* [ndim] file dimensions */, 
		      int argc, char**argv /* command-line arguments */)
/*< split the input file along the specified axis
  and generate parallel system commands >*/
{
    char **commands, **splitkey, okey[5], dkey[5];
    char *arg, *eq, *cmdline, *iname=NULL, *oname=NULL, *splitname;
    off_t i2, *splitsize1=NULL, *splitsize2=NULL, left, nbuf;
    float d, o, di, oi;
    int job, jobs, bigjobs, w, split, i, j, len, chunk, skip;
    sf_file *splitfile=NULL, in=NULL;
    FILE *ifile=NULL, *ofile=NULL;

    if (axis > ndim) axis=ndim;

    snprintf(nkey,5,"n%d",axis);
    snprintf(okey,5,"o%d",axis);
    snprintf(dkey,5,"d%d",axis);
    if (!sf_histfloat(inp,okey,&o)) o=0.;
    if (!sf_histfloat(inp,dkey,&d)) d=1.;
    
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
    
    splitkey = (char**) sf_alloc(argc,sizeof(char*));
    splitargc = 0;

    j=0;
    for (i=1; i < argc; i++) {
	arg = argv[i];
	if (strncmp(arg,"input=",6) &&
	    strncmp(arg,"output=",7) &&
	    strncmp(arg,"split=",6) &&
	    strncmp(arg,"join=",5)) {
	    if ('_'==arg[0]) {
		eq  = strchr(arg,'=');
		if (NULL == eq) sf_error("Wrong parameter \"%s\"",arg);

		len = eq-arg;
		splitkey[splitargc] = sf_charalloc(len+1);
		strncpy(splitkey[splitargc],arg,len);
		splitkey[splitargc][len]='\0';

		splitargc++;
	    } else {
		len = strlen(arg);
		if (j+len > SF_CMDLEN-2) sf_error("command line is too long");

		strncpy(command+j,arg,len);
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
	    splitfile[i] = sf_input(splitkey[i]);
	    ndim = sf_largefiledims (splitfile[i],n);
	    sizes(splitfile[i],axis,ndim,n,splitsize1+i,splitsize2+i);
	    
	    if (n[axis] != split) 
		sf_error("Wrong dimensions in file %s",splitkey[i]);
	}
	    
	spnames = (char***) sf_alloc(jobs,sizeof(char**));
	for (job=0; job < jobs; job++) {
	    spnames[job] = (char**) sf_alloc(splitargc,sizeof(char*));
	}
    }

    inames = (char**) sf_alloc(jobs,sizeof(char*));
    onames = (char**) sf_alloc(jobs,sizeof(char*));
    commands =  (char**) sf_alloc(jobs,sizeof(char*));

    for (job=0; job < jobs; job++) {
	commands[job] = sf_charalloc(SF_CMDLEN);
	cmdline = commands[job];

	if (job < bigjobs) {
	    chunk = w;
	    skip = job*w;
	} else {
	    chunk = w-1;
	    skip = bigjobs*w+(job-bigjobs)*chunk;
	}

	ifile = sf_tempfile(&iname,"w+b");
	ofile = sf_tempfile(&oname,"w+b");
	fclose(ofile);

	inames[job]=iname;
	onames[job]=oname;

	in = sf_output(iname);
	fclose(ifile);

	sf_putint(in,nkey,chunk);
	sf_putfloat(in,okey,o+skip*d);
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
	    ifile = sf_tempfile(&splitname,"w+b");

	    spnames[job][i] = splitname;

	    in = sf_output(splitname);
	    fclose(ifile);

	    if (!sf_histfloat(splitfile[i],okey,&oi)) oi=0.;
	    if (!sf_histfloat(splitfile[i],dkey,&di)) di=1.;

	    sf_putint(in,nkey,chunk);
	    sf_putfloat(in,okey,oi+skip*di);
	    sf_fileflush(in,splitfile[i]);
	    sf_setform(in,SF_NATIVE);

	    for (i2=0; i2 < splitsize2[i]; i2++) {
		sf_seek(splitfile[i],splitsize1[i]*(i2*split+skip),SEEK_SET);

		nbuf=BUFSIZ;
		for (left=chunk*splitsize1[i]; left > 0; left -= nbuf) {
		    if (nbuf > left) nbuf=left;
			
		    sf_charread(buffer,nbuf,splitfile[i]);
		    sf_charwrite(buffer,nbuf,in);
		}
	    }

	    sf_fileclose(in);

	    snprintf(cmdline,SF_CMDLEN,"%s %s=%s",
		     splitcommand,splitkey[i]+1,splitname);
	    strncpy(splitcommand,cmdline,SF_CMDLEN);
	}	

	snprintf(cmdline,SF_CMDLEN,"%s %s < %s > %s",
		 command,splitcommand,iname,oname);
    }

    return commands;
}

void sf_out(sf_file out          /* output file */, 
		  int axis             /* join axis */,
		  const char *iname    /* name of the input file */)
/*< prepare output >*/
{
    char *oname, cmdline[SF_CMDLEN];
    int ndim;
    off_t n[SF_MAX_DIM];
    sf_file inp;
    FILE *ofile;

    ofile = sf_tempfile(&oname,"w+b");
    
    snprintf(cmdline,SF_CMDLEN,"%s %s --dryrun=y < %s > %s",
	     command,splitcommand,iname,oname);
    sf_system(cmdline);
    
    inp = sf_input(oname);
    ndim = sf_largefiledims (inp,n);
    
    if (axis > ndim) axis=ndim;

    snprintf(nkey,5,"n%d",axis);
    axis--;
    sizes(inp,axis,ndim,n,&size1,&size2);
    
    sf_setformat(out,sf_histstring(inp,"data_format"));
    sf_fileflush(out,inp);
    
    sf_setform(out,SF_NATIVE);
    sf_rm(oname,true,false,false);
}

void sf_join(sf_file out, int job)
/*< join outputs >*/
{
    int i, chunk;
    off_t i2, left, nbuf;
    char *iname, *oname;
    sf_file in;

    iname = inames[job];
    oname = onames[job];

    in = sf_input(oname);
    sf_histint(in,nkey,&chunk);
    sf_setform(in,SF_NATIVE);
	    
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

    for (i=0; i < splitargc; i++) {
	iname = spnames[job][i];
	sf_rm(iname,true,false,false);
    }
}
