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
static char ***inpnames=NULL, buffer[BUFSIZ], nkey[5];
static off_t size1, size2;
static int inpargc, outargc;
static sf_file *ins;

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
		int *tasks           /* number of tasks */,
		int ndim, off_t *n   /* [ndim] file dimensions */, 
		int argc, char**argv /* command-line arguments */)
/*< split the input file along the specified axis
  and generate parallel system commands >*/
{
    char **commands, **splitinp, **splitout, okey[5], dkey[5], splitcmd[SF_CMDLEN];
    char *arg, *eq, *cmdline, *iname=NULL, *oname=NULL, *splitname;
    off_t i2, *splitsize1=NULL, *splitsize2=NULL, left, nbuf;
    float d, o, di, oi;
    int job, jobs, bigjobs, w, split, i, j, k, len, chunk, skip;
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
    if (jobs < split) {
	w = (int) (1+((float) split)/jobs);
    } else {
	w = 1;
	jobs = split;
    }
    bigjobs = split - jobs*(w-1);
    *tasks = jobs;
    
    splitinp = (char**) sf_alloc(argc,sizeof(char*));
    splitout = (char**) sf_alloc(argc,sizeof(char*));
    outargc = 0;
    inpargc = 0;

    k=j=0;
    for (i=1; i < argc; i++) {
	arg = argv[i];
	if (strncmp(arg,"--input=",8) &&
	    strncmp(arg,"--output=",9) &&
	    strncmp(arg,"split=",6) &&
	    strncmp(arg,"join=",5)) {

	    len = strlen(arg);
	    if (j+len > SF_CMDLEN-2) sf_error("command line is too long");

	    if (len > 2 && '_'==arg[0] && '_'==arg[1]) {
		/* starting a parameter with two underscores signifies an output
		 * file to split */

		eq  = strchr(arg,'=');
		if (NULL == eq) sf_error("Wrong parameter \"%s\"",arg);

		len = eq-arg;
		splitout[outargc] = sf_charalloc(len+1);
		strncpy(splitout[outargc],arg,len);
		splitout[outargc][len]='\0';

		outargc++;

		len = strlen(arg)-2;
		if (k+len > SF_CMDLEN-2) sf_error("command line is too long");
		strncpy(splitcommand+k,arg+2,len);
		splitcommand[k+len]=' ';
		k += len+1;
	    } else if ('_'==arg[0]) {
		/* starting a parameter with one underscore signifies an input
		 * file to split */

		eq  = strchr(arg,'=');
		if (NULL == eq) sf_error("Wrong parameter \"%s\"",arg);

		len = eq-arg;
		splitinp[inpargc] = sf_charalloc(len+1);
		strncpy(splitinp[inpargc],arg,len);
		splitinp[inpargc][len]='\0';

		inpargc++;

		len = strlen(arg)-1;
		if (k+len > SF_CMDLEN-2) sf_error("command line is too long");
		strncpy(splitcommand+k,arg+1,len);
		splitcommand[k+len]=' ';
		k += len+1;
	    } else {
		strncpy(command+j,arg,len);
		command[j+len]=' ';
		j += len+1;
	    }
	}
    }
    command[j]='\0';
    splitcommand[k]='\0';

    if (0 < inpargc) { /* files to split other than input */
	splitsize1 = sf_largeintalloc(inpargc);
	splitsize2 = sf_largeintalloc(inpargc);
	splitfile = (sf_file*) sf_alloc(inpargc,sizeof(sf_file));

	for (i=0; i < inpargc; i++) {
	    splitfile[i] = sf_input(splitinp[i]);
	    ndim = sf_largefiledims (splitfile[i],n);
	    sizes(splitfile[i],axis,ndim,n,splitsize1+i,splitsize2+i);
	    
	    if (n[axis] != split) 
		sf_error("Wrong dimensions in file %s",splitinp[i]);
	}
	    
	inpnames = (char***) sf_alloc(jobs,sizeof(char**));
	for (job=0; job < jobs; job++) {
	    inpnames[job] = (char**) sf_alloc(inpargc,sizeof(char*));
	}
    }

    inames = (char**) sf_alloc(jobs,sizeof(char*));
    onames = (char**) sf_alloc(jobs,sizeof(char*));
    commands =  (char**) sf_alloc(jobs,sizeof(char*));

    for (job=0; job < jobs; job++) {
	commands[job] = sf_charalloc(jobs*SF_CMDLEN+1);
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
	sf_settype(in,sf_gettype(inp));
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

	splitcmd[0]='\0';
	for (i=0; i < inpargc; i++) {
	    ifile = sf_tempfile(&splitname,"w+b");

	    inpnames[job][i] = splitname;

	    in = sf_output(splitname);
	    fclose(ifile);

	    if (!sf_histfloat(splitfile[i],okey,&oi)) oi=0.;
	    if (!sf_histfloat(splitfile[i],dkey,&di)) di=1.;

	    sf_putint(in,nkey,chunk);
	    sf_putfloat(in,okey,oi+skip*di);
	    sf_settype(in,sf_gettype(splitfile[i]));
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

	    snprintf(cmdline,jobs*SF_CMDLEN,"%s %s=%s",splitcmd,splitinp[i]+1,splitname);
	    strncpy(splitcmd,cmdline,SF_CMDLEN);
	}	

	snprintf(cmdline,jobs*SF_CMDLEN,"%s %s < %s > %s",command,splitcmd,iname,oname);
    }

    return commands;
}

void sf_out(sf_file out        /* output file */,
	    int jobs           /* number of jobs */,
	    int axis           /* join axis */,
	    const char *iname  /* name of the input file */)
/*< prepare output >*/
{
    char *oname, *cmdline;
    int ndim, job;
    off_t n[SF_MAX_DIM];
    sf_file inp;
    FILE *ofile=NULL;

    ofile = sf_tempfile(&oname,"w+b");
    fclose(ofile);

    cmdline = sf_charalloc(jobs*SF_CMDLEN+1);
    snprintf(cmdline,jobs*SF_CMDLEN,"%s %s --dryrun=y < %s > %s",
	     command,splitcommand,iname,oname);
    sf_system(cmdline);
    free(cmdline);
    
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

    ins = (sf_file*) sf_alloc(jobs,sizeof(sf_file));

    for (job=0; job < jobs; job++) {
	oname = onames[job];
	ins[job] = sf_input(oname);
    }
}

static void cleanup(int jobs)
/* remove temporary files */
{
    int i, job;
    char *oname;

    for (job=0; job < jobs; job++) {
	sf_fileclose(ins[job]);
	sf_rm(inames[job],true,false,false);
	sf_rm(onames[job],true,false,false);
	for (i=0; i < inpargc; i++) {
	    oname = inpnames[job][i];
	    sf_rm(oname,true,false,false);
	}
    }

    free(ins);
}

void sf_join(sf_file out /* output file */,
	     int axis    /* axis to join */,
	     int jobs    /* number of jobs */)
/*< join outputs >*/
{
    int i, job, dim, esize, *naxis, ni;    
    off_t i2, left, nbuf, n[SF_MAX_DIM], n1, n2;
    char key[12], buf[BUFSIZ];
    sf_file in;

    dim = sf_largefiledims(ins[0],n);
    esize = sf_esize(ins[0]);
    naxis = sf_intalloc(jobs);
    
    n1=1;
    n2=1;
    for (i=1; i <= dim; i++) {
	if      (i < axis) n1 *= n[i-1];
	else if (i > axis) n2 *= n[i-1];
    }

    naxis[0] = n[axis-1];
    for (job=1; job < jobs; job++) {
	in = ins[job];
	
	for (i=1; i <= dim; i++) {
	    (void) snprintf(key,3,"n%d",i%10u);
	    if (!sf_histint(in,key,&ni) || (i != axis && ni != n[i-1]))
		sf_error("%s mismatch: need %d",key,(int) n[i-1]);
	    if (i == axis) naxis[job] = ni;
	}
	if (axis > dim) {
	    (void) snprintf(key,3,"n%d",axis);
	    if (!sf_histint(in,key,naxis+job)) naxis[job]=1;
	}
    } 
    
    for (i2=0; i2 < n2; i2++) {
	for (job=0; job < jobs; job++) {
	    
	    in = ins[job];
	    sf_setform(in,SF_NATIVE);
	    
	    for (left = n1*naxis[job]*esize; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ: left;
		sf_charread (buf,nbuf,in);
		sf_charwrite (buf,nbuf,out);
	    }
	}
    }

    free(naxis);
    cleanup(jobs);
}

void sf_add(sf_file out /* output file */,
	    int jobs    /* number of jobs */)
/*< add outputs >*/
{
    int *ibuf=NULL, job, i;
    float *fbuf=NULL;
    char buffout[BUFSIZ];
    sf_complex *cbuf=NULL;
    sf_datatype type;
    size_t nbuf=BUFSIZ;
    off_t nsiz;

    type = sf_gettype(out);

    switch(type) {
	case SF_FLOAT:
	    fbuf = (float*) buffout;
	    nbuf /= sizeof(float);
	    for (i=0; i < nbuf; i++) {
		fbuf[i] = 0.0f;
	    }
	    break;
	default:
	    sf_error("wrong type");
	    break;
    }

    for (nsiz = size2; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;

	for (job=0; job < jobs; job++) {
	    switch(type) {
		case SF_FLOAT:
		    sf_floatread((float*) buffer,nbuf,ins[job]);
		    for (i=0; i < nbuf; i++) {
			if (job) {
			    fbuf[i] += ((float*) buffer)[i];
			} else {
			    fbuf[i] = ((float*) buffer)[i];
			}
		    }
		    break;
		default:
		    sf_error("wrong type");
		    break;  
	    }
	}

	switch(type) {
	    case SF_FLOAT:
		sf_floatwrite(fbuf,nbuf,out);
		break;
	    case SF_COMPLEX:
		sf_complexwrite(cbuf,nbuf,out);
		break;
	    case SF_INT:
		sf_intwrite(ibuf,nbuf,out);
		break;
	    default:
		sf_error("wrong type");
		break;
	}
    }

    cleanup(jobs);
}

