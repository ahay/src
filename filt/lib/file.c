/* Main operations with RSF files. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include <stdio.h>
/*^*/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __GNUC__
#ifndef alloca
#define alloca __builtin_alloca
#endif
#else /* not GNU C  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || defined(hpux) || defined(__hpux)
#include <alloca.h>
#endif
#endif

#include <limits.h>

#include <sys/stat.h>
#include <sys/param.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>

#include <rpc/types.h>
/* #include <rpc/rpc.h> */
#include <rpc/xdr.h>

#include "_defs.h"
#include "file.h"
#include "getpar.h"
#include "alloc.h"
#include "error.h"
#include "simtab.h"
#include "komplex.h"

#include "_bool.h"
#include "c99.h"
/*^*/

/* BSD - MAXNAMELEN, Posix - NAME_MAX */
#ifndef NAME_MAX
#ifdef MAXNAMELEN
#define	NAME_MAX MAXNAMELEN
#else
#ifdef FILENAME_MAX
#define NAME_MAX FILENAME_MAX
#endif
#endif
#endif

#ifndef _sf_file_h

#define SF_MAX_DIM 9
/*^*/

typedef struct sf_File *sf_file;
/*^*/

typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;
/*^*/

#endif

struct sf_File {
    FILE *stream, *head; 
    char *dataname, *buf;
    sf_simtab pars;
    XDR xdr;
    enum xdr_op op;
    sf_datatype type;
    sf_dataform form;
    bool pipe, rw, dryrun;
};

/*@null@*/ static sf_file infile=NULL;
static const int tabsize=10;
/*@null@*/ static char *aformat = NULL;
static size_t aline=8;

static bool getfilename (FILE *fd, char *filename);
static char* getdatapath (void);
static char* gettmpdatapath (void);
static bool readpathfile (const char* filename, char* datapath);

sf_file sf_input (/*@null@*/ const char* tag)
/*< Create an input file structure >*/
{
    int esize;
    sf_file file;
    char *filename, *format;
    size_t len;
    extern off_t ftello (FILE *stream);

    file = (sf_file) sf_alloc(1,sizeof(*file));
    
    if (NULL == tag || 0 == strcmp(tag,"in")) {
	file->stream = stdin;
	filename = NULL;
    } else {
	filename = sf_getstring (tag);
	if (NULL == filename) {
	    len = strlen(tag)+1;
	    filename = sf_charalloc(len);
	    memcpy(filename,tag,len);
	}

	file->stream = fopen(filename,"r");
	if (NULL == file->stream) 
	    sf_error ("%s: Cannot read header file %s:",__FILE__,filename);
    }
/*    setbuf(file->stream,file->buf); */

    /* create a parameter table */
    file->pars = sf_simtab_init (tabsize);
    file->head = tmpfile();

    /* read the parameter table from the file */
    sf_simtab_input (file->pars,file->stream,file->head);

    if (NULL == filename) {
	infile = file;
    } else {
	free (filename);
    }

    filename = sf_histstring(file,"in");
    if (NULL == filename) sf_error ("%s: No in= in file %s",__FILE__,tag);
    len = strlen(filename)+1;
    file->dataname = sf_charalloc(len);
    memcpy(file->dataname,filename,len);

    if (0 != strcmp(filename,"stdin")) {
	file->stream = freopen(filename,"rb",file->stream);
	if (NULL == file->stream) 
	    sf_error("%s: Cannot read data file %s:",__FILE__,filename);
    }
    free (filename);

    file->pipe = (bool) (-1 == ftello(file->stream));
    if (file->pipe && ESPIPE != errno) 
	sf_error ("%s: pipe problem:",__FILE__);

    file->op = XDR_DECODE;
    file->buf = NULL;

    format = sf_histstring(file,"data_format");
    if (NULL == format) {
	if (!sf_histint(file,"esize",&esize) || 0 != esize)
	    sf_error ("%s: Unknown format in %s",__FILE__,tag);
	sf_setformat(file,"ascii_float");
    } else {    
	sf_setformat(file,format);
	free (format);
    }

    return file;
}

sf_file sf_output (/*@null@*/ const char* tag)
/*< Create an output file structure.
---
Should do output after sf_input. >*/
{
    sf_file file;
    char *headname, *dataname, *path, *name, *format;
    size_t namelen;
    extern int mkstemp (char *template);
    extern off_t ftello (FILE *stream);

    file = (sf_file) sf_alloc(1,sizeof(*file));

    if (NULL == tag || 0 == strcmp(tag,"out")) {
	file->stream = stdout;
	headname = NULL;
    } else {
	headname = sf_getstring (tag);
	if (NULL == headname) {
	    namelen = strlen(tag)+1;
	    headname = sf_charalloc (namelen);
	    memcpy(headname,tag,namelen);
	}

	file->stream = fopen(headname,"w");
	if (NULL == file->stream) 
	    sf_error ("%s: Cannot write to header file %s:",__FILE__,headname);
    }
/*    setbuf(file->stream,file->buf); */

    file->pars = sf_simtab_init (tabsize);
    file->head = NULL;

    file->pipe = (bool) (-1 == ftello(file->stream));
    if (file->pipe && ESPIPE != errno) 
	sf_error ("%s: pipe problem:",__FILE__);

    dataname = sf_getstring("out");

    if (file->pipe) {
	file->dataname = sf_charalloc (7);
	memcpy(file->dataname,"stdout",7);
    } else if (NULL == dataname) {
	path = getdatapath();
	file->dataname = sf_charalloc (PATH_MAX+NAME_MAX+1);
	strcpy (file->dataname,path);
	name = file->dataname+strlen(path);
	free (path);
	if (getfilename (file->stream,name)) {
	    namelen = strlen(file->dataname);
	    file->dataname[namelen]='@';
	    file->dataname[namelen+1]='\0';
	} else { /* invent a name */
	    sprintf(name,"%sXXXXXX",sf_getprog());
	    (void) close(mkstemp(file->dataname));
	    /* (void) unlink(file->dataname); */
	    /* old code for named pipes below */
	    /*
	    if (NULL == headname &&
		-1L == fseek(file->stream,0L,SEEK_CUR) &&
		ESPIPE == errno && 
		0 != mkfifo (file->dataname, S_IRUSR | S_IWUSR))
		sf_error ("%s: Cannot make a pipe %s:",
			  __FILE__,file->dataname);
	    */
	}  
    } else {
	namelen = strlen(dataname)+1;
	file->dataname = sf_charalloc (namelen);
	memcpy(file->dataname,dataname,namelen);
	free (dataname);
    }

    sf_putstring(file,"in",file->dataname);    

    file->op = XDR_ENCODE;
    file->buf = NULL;

    if (NULL != infile && 
	NULL != (format = sf_histstring(infile,"data_format"))) {
	sf_setformat(file,format);
	free (format);
    } else {
	sf_setformat(file,"native_float");
    }

    if (NULL != headname) free(headname);

    if (!sf_getbool("readwrite",&(file->rw))) file->rw=false;
    if (!sf_getbool("dryrun",&(file->dryrun))) file->dryrun=false;

    return file;
}

sf_datatype sf_gettype (sf_file file)
/*< return file type >*/
{
    return file->type;
}

sf_dataform sf_getform (sf_file file)
/*< return file form >*/
{
    return file->form;
}

size_t sf_esize(sf_file file)
/*< return element size >*/
{
    sf_datatype type;

    type = sf_gettype (file);
    switch (type) {
	case SF_FLOAT: 
	    return sizeof(float);
	    break;
	case SF_INT:
	    return sizeof(int);
	    break;
	case SF_COMPLEX:
	    return 2*sizeof(float);
	    break;
	default:
	    return sizeof(char);
	    break;
    }
}


void sf_settype (sf_file file, sf_datatype type)
/*< set file type >*/
{
    file->type = type;
    sf_putint(file,"esize",(int) sf_esize(file));
}

void sf_setpars (sf_file file)
/*< change parameters to those from the command line >*/
{
    char *in;

    in = sf_histstring(file,"in");
    sf_simtab_close(file->pars);
    file->pars = sf_getpars();

    if (NULL != in) {
	sf_putstring(file,"in",in);
	free(in);
    }
}

size_t sf_bufsiz(sf_file file)
/*< return buffer size for efficient I/O >*/
{
    size_t bufsiz;
    int filedes;
    struct stat buf;

    filedes = fileno(file->stream);
    if (fstat(filedes,&buf)) {
	bufsiz = BUFSIZ;
    } else {
	bufsiz = buf.st_blksize;
	if (bufsiz <=0) bufsiz = BUFSIZ;
    } 
    
    return bufsiz;
}	

void sf_setform (sf_file file, sf_dataform form)
/*< set file form >*/
{
    size_t bufsiz;

    file->form = form;
    
    switch(form) {
	case SF_ASCII:
	    if (NULL != file->buf) {
		free (file->buf);
		file->buf = NULL;
	    }
	    sf_putint(file,"esize",0); /* for compatibility with SEPlib */
	    break;
	case SF_XDR:
	    if (NULL == file->buf) {
		bufsiz = sf_bufsiz(file);
		file->buf = sf_charalloc(bufsiz);
		xdrmem_create(&(file->xdr),file->buf,bufsiz,file->op);
	    }
	    break;
	case SF_NATIVE:
	default:
	    if (NULL != file->buf) {
		free (file->buf);
		file->buf = NULL;
	    }
	    break;
    }
}

void sf_setformat (sf_file file, const char* format)
/*< Set file format.
---
format has a form "form_type", i.e. native_float, ascii_int, etc.
>*/
{
    if (NULL != strstr(format,"float")) {
	sf_settype(file,SF_FLOAT);
    } else if (NULL != strstr(format,"int")) {
	sf_settype(file,SF_INT);
    } else if (NULL != strstr(format,"complex")) {
	sf_settype(file,SF_COMPLEX);
    } else if (NULL != strstr(format,"uchar") || 
	       NULL != strstr(format,"byte")) {
	sf_settype(file,SF_UCHAR);	
    } else {
	sf_settype(file,SF_CHAR);
    }
    
    if (0 == strncmp(format,"ascii_",6)) {
	sf_setform(file,SF_ASCII);
    } else if (0 == strncmp(format,"xdr_",4)) {
	sf_setform(file,SF_XDR);
    } else {
	sf_setform(file,SF_NATIVE);
    }
}
    
static bool getfilename (FILE* fp, char *filename)
/* Finds filename of an open file from the file descriptor.

Unix-specific and probably non-portable. */
{
    DIR* dir;
    struct stat buf;
    struct dirent *dirp;
    bool success;

    dir = opendir(".");
    if (NULL == dir) sf_error ("%s: cannot open current directory:",__FILE__);
    
    if(0 > fstat(fileno(fp),&buf)) sf_error ("%s: cannot run fstat:",__FILE__);
    success = false;
    
    while (NULL != (dirp = readdir(dir))) {
	if (dirp->d_ino == buf.st_ino) { /* non-portable */
	    strcpy(filename,dirp->d_name);
	    success = true;
	    break;
	}
    }

    closedir(dir);

    return success;
}

static char* gettmpdatapath (void) 
/* Finds temporary datapath.

Datapath rules:
1. check tmpdatapath= on the command line
2. check TMPDATAPATH environmental variable
3. check .tmpdatapath file in the current directory
4. check .tmpdatapath in the home directory
5. return NULL
*/
{
    char *path, *penv, *home, file[PATH_MAX];
    
    path = sf_getstring ("tmpdatapath");
    if (NULL != path) return path;

    path = sf_charalloc(PATH_MAX+1);

    penv = getenv("TMPDATAPATH");
    if (NULL != penv) {
	strncpy(path,penv,PATH_MAX);
	return path;
    }

    if (readpathfile (".tmpdatapath",path)) return path;
    
    home = getenv("HOME");
    if (NULL != home) {
	(void) snprintf(file,PATH_MAX,"%s/.tmpdatapath",home);
	if (readpathfile (file,path)) return path;
    }

    return NULL;
}

static char* getdatapath (void) 
/* Finds datapath.

Datapath rules:
1. check datapath= on the command line
2. check DATAPATH environmental variable
3. check .datapath file in the current directory
4. check .datapath in the home directory
5. use '.' (not a SEPlib behavior)
*/
{
    char *path, *penv, *home, file[PATH_MAX];
    
    path = sf_getstring ("datapath");
    if (NULL != path) return path;

    path = sf_charalloc(PATH_MAX+1);

    penv = getenv("DATAPATH");
    if (NULL != penv) {
	strncpy(path,penv,PATH_MAX);
	return path;
    }

    if (readpathfile (".datapath",path)) return path;
    
    home = getenv("HOME");
    if (NULL != home) {
	(void) snprintf(file,PATH_MAX,"%s/.datapath",home);
	if (readpathfile (file,path)) return path;
    }

    path = sf_charalloc(3);
    strncpy(path,"./",3);

    return path;
}

static bool readpathfile (const char* filename, char* datapath) 
/* find datapath from the datapath file */
{
    FILE *fp;
    char host[PATH_MAX], *thishost, path[PATH_MAX];

    fp = fopen(filename,"r");
    if (NULL == fp) return false;

    if (0 >= fscanf(fp,"datapath=%s",datapath))
	sf_error ("No datapath found in file %s",filename);

    thishost = sf_gethost();

    while (0 < fscanf(fp,"%s datapath=%s",host,path)) {
	if (0 == strcmp(host,thishost)) {
	    strcpy(datapath,path);
	    break;
	}
    }

    (void) fclose (fp);
    return true;
}

void sf_fileclose (sf_file file) 
/*< close a file and free allocated space >*/
{
    if (file->stream != stdin && 
	file->stream != stdout && 
	file->stream != NULL) (void) fclose (file->stream);
    if (NULL != file->pars) sf_simtab_close (file->pars);
    if (NULL != file->buf) {
	free (file->buf);
	file->buf = NULL;
    }
    if (NULL != file->dataname) free (file->dataname);
    free (file);
}

bool sf_histint (sf_file file, const char* key,/*@out@*/ int* par)
/*< read an int parameter from file >*/ 
{
    return sf_simtab_getint (file->pars,key,par);
}

bool sf_histints (sf_file file, const char* key,/*@out@*/ int* par,size_t n) 
/*< read an int array of size n parameter from file >*/ 
{
    return sf_simtab_getints (file->pars,key,par, n);
}

bool sf_histlongint (sf_file file, const char* key,/*@out@*/ long int* par)
/*< read a long int parameter from file >*/ 
{
    return sf_simtab_getlongint (file->pars,key,par);
}

bool sf_histfloat (sf_file file, const char* key,/*@out@*/ float* par) 
/*< read a float parameter from file >*/
{
    return sf_simtab_getfloat (file->pars,key,par);
}

bool sf_histdouble (sf_file file, const char* key,/*@out@*/ double* par) 
/*< read a float parameter from file >*/
{
    return sf_simtab_getdouble (file->pars,key,par);
}

bool sf_histfloats (sf_file file, const char* key,
		    /*@out@*/ float* par,size_t n) 
/*< read a float array of size n parameter from file >*/ 
{
    return sf_simtab_getfloats (file->pars,key,par, n);
}

bool sf_histbool (sf_file file, const char* key,/*@out@*/ bool* par) 
/*< read a bool parameter from file >*/
{
    return sf_simtab_getbool (file->pars,key,par);
}

bool sf_histbools (sf_file file, const char* key,
		   /*@out@*/ bool* par, size_t n) 
/*< read a bool array of size n parameter from file >*/ 
{
    return sf_simtab_getbools (file->pars,key,par, n);
}

char* sf_histstring (sf_file file, const char* key) 
/*< read a string parameter from file (returns NULL on failure) >*/ 
{
    return sf_simtab_getstring (file->pars,key);
}

void sf_fileflush (sf_file file, sf_file src)
/*< outputs parameter to a file (initially from source src)
---
Prepares file for writing binary data >*/ 
{
    time_t tm;
    char line[BUFSIZ];
 
    /* if already flushed, do nothing */
    if (NULL == file->dataname) return;
 
    if (NULL != src && NULL != src->head) {
	rewind(src->head);

	while (NULL != fgets(line,BUFSIZ,src->head)) {
	    fputs(line,file->stream);
	}
    }
    
    tm = time (NULL);
    if (0 >= fprintf(file->stream,"%s\t%s:\t%s@%s\t%s\n",
		     sf_getprog(),
		     sf_getcdir(),
		     sf_getuser(),
		     sf_gethost(),
		     ctime(&tm)))
	sf_error ("%s: cannot flush header:",__FILE__);

    switch (file->type) {
	case SF_FLOAT: 
	    switch (file->form) {
		case SF_ASCII:		    
		    sf_putstring(file,"data_format","ascii_float");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_float");
		    break;
		default:
		    sf_putstring(file,"data_format","native_float");
		    break;
	    }
	    break;
	case SF_COMPLEX:
	    switch (file->form) {
		case SF_ASCII:		    
		    sf_putstring(file,"data_format","ascii_complex");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_complex");
		    break;
		default:
		    sf_putstring(file,"data_format","native_complex");
		    break;
	    }
	    break;
	case SF_INT:
	    switch (file->form) {
		case SF_ASCII:		    
		    sf_putstring(file,"data_format","ascii_int");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_int");
		    break;
		default:
		    sf_putstring(file,"data_format","native_int");
		    break;
	    }
	    break;
	case SF_UCHAR:
	    switch (file->form) {
		case SF_ASCII:		    
		    sf_putstring(file,"data_format","ascii_uchar");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_uchar");
		    break;
		default:
		    sf_putstring(file,"data_format","native_uchar");
		    break;
	    }
	    break;
	case SF_CHAR:
	    switch (file->form) {
		case SF_ASCII:		    
		    sf_putstring(file,"data_format","ascii_char");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_char");
		    break;
		default:
		    sf_putstring(file,"data_format","native_char");
		    break;
	    }
	    break;
	default:
	    sf_putstring(file,"data_format",
			 (NULL==file->buf)? "native_byte":"xdr_byte");
	    break;
    }    
    
    sf_simtab_output(file->pars,file->stream);

    if (0==strcmp(file->dataname,"stdout")) { 
	/* keep stream, write the header end code */
	fprintf(file->stream,"\tin=\"stdin\"\n\n%c%c%c",
		SF_EOL,SF_EOL,SF_EOT);
    } else {
	file->stream = freopen(file->dataname,file->rw? "w+b":"wb",file->stream);       
	if (NULL == file->stream) 
	    sf_error ("%s: Cannot write to data file %s:",
		      __FILE__,file->dataname);	
    }
    
    free (file->dataname);
    file->dataname=NULL;

    if (file->dryrun) exit(0);
}

void sf_putint (sf_file file, const char* key, int par)
/*< put an int parameter to a file >*/
{
    char val[256];

    snprintf(val,256,"%d",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putints (sf_file file, const char* key, const int* par, size_t n)
/*< put an int array of size n parameter to a file >*/
{
    int i;
    char val[1024], *v;

    v = val;
    for (i=0; i < n-1; i++) {
	v += snprintf(v,1024,"%d,",par[i]);
    }
    snprintf(v,1024,"%d",par[n-1]);

    sf_simtab_enter (file->pars,key,val);
}

void sf_putlongint (sf_file file, const char* key, long int par)
/*< put a long int parameter to a file >*/
{
    char val[256];

    snprintf(val,256,"%ld",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putfloat (sf_file file, const char* key,float par)
/*< put a float parameter to a file >*/
{
    char val[256];

    snprintf(val,256,"%g",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putstring (sf_file file, const char* key,const char* par)
/*< put a string parameter to a file >*/
{
    char *val;
    
    val = (char*) alloca(strlen(par)+3); 
    sprintf(val,"\"%s\"",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putline (sf_file file, const char* line)
/*< put a string line to a file >*/
{
    if (0 >= fprintf(file->stream,"\t%s\n",line))
	sf_error ("%s: cannot put line '%s':",__FILE__,line);
}

void sf_setaformat (const char* format /* number format (.i.e "%5g") */, 
		    int line /* numbers in line */ )
/*< Set format for ascii output >*/
{
    size_t len;

    if (NULL != aformat) free (aformat);
    if (NULL != format) {
	len = strlen(format)+1;
	aformat = sf_charalloc(len);
	memcpy(aformat,format,len);
    } else {
	aformat = NULL;
    }
    aline = (size_t) line;
}

void sf_complexwrite (sf_complex* arr, size_t size, sf_file file)
/*< write a complex array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;
    sf_complex c;

    if (NULL != file->dataname) sf_fileflush (file,infile);
    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    c = arr[i];
		    if (EOF==fprintf(file->stream,
				     (NULL != aformat)? 
				     aformat:"%g %gi ",
				     crealf(c),cimagf(c)))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(sf_complex);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(sf_complex),size,file->stream);
	    break;
    }
}

void sf_complexread (/*@out@*/ sf_complex* arr, size_t size, sf_file file)
/*< read a complex array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;
    float re, im;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		if (EOF==fscanf(file->stream,"%g %gi",&re,&im))
		    sf_error ("%s: trouble reading ascii:",__FILE__);
		arr[i]=sf_cmplx(re,im);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(sf_complex);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(sf_complex),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

void sf_charwrite (char* arr, size_t size, sf_file file)
/*< write a char array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;

    if (NULL != file->dataname) sf_fileflush (file,infile);

    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    if (EOF==fputc(arr[i],file->stream))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    buf = arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_opaque(&(file->xdr),buf-left,nbuf))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(char),size,file->stream);
	    break;
    }
}

void sf_ucharwrite (unsigned char* arr, size_t size, sf_file file)
/*< write an unsigned char array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;

    if (NULL != file->dataname) sf_fileflush (file,infile);

    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    if (EOF==fputc(arr[i],file->stream))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_opaque(&(file->xdr),buf-left,nbuf))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(unsigned char),size,file->stream);
	    break;
    }
}

void sf_charread (/*@out@*/ char* arr, size_t size, sf_file file)
/*< read a char array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;
    int c;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		c=fgetc(file->stream);
		if (EOF==c)
		    sf_error ("%s: trouble reading ascii:",__FILE__);
		arr[i]= (char) c;
	    }
	    break;
	case SF_XDR:
	    buf = arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_opaque(&(file->xdr),buf-left,nbuf))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(char),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

void sf_ucharread (/*@out@*/ unsigned char* arr, size_t size, sf_file file)
/*< read a uchar array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;
    int c;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		c=fgetc(file->stream);
		if (EOF==c)
		    sf_error ("%s: trouble reading ascii:",__FILE__);
		arr[i]= (unsigned char) c;
	    }
	    break;
	case SF_XDR:
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_opaque(&(file->xdr),buf-left,nbuf))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(unsigned char),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

void sf_intwrite (int* arr, size_t size, sf_file file)
/*< write an int array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;

    if (NULL != file->dataname) sf_fileflush (file,infile);
    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    if (EOF==fprintf(file->stream,
				     (NULL != aformat)? aformat:"%d ",
				     arr[i]))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(int);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(int),sizeof(int),
				(xdrproc_t) xdr_int))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(int),size,file->stream);
	    break;
    }
}

void sf_intread (/*@out@*/ int* arr, size_t size, sf_file file)
/*< read an int array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		if (EOF==fscanf(file->stream,"%d",arr+i))
		    sf_error ("%s: trouble reading ascii:",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(int);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(int),sizeof(int),
				(xdrproc_t) xdr_int))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(int),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

void sf_floatwrite (float* arr, size_t size, sf_file file)
/*< write a float array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;

    if (NULL != file->dataname) sf_fileflush (file,infile);
    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    if (EOF==fprintf(file->stream,
				     (NULL != aformat)? aformat:"%g ",
				     arr[i]))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(float);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(float),size,file->stream);
	    break;
    }
}

void sf_floatread (/*@out@*/ float* arr, size_t size, sf_file file)
/*< read a float array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		if (EOF==fscanf(file->stream,"%g",arr+i))
		    sf_error ("%s: trouble reading ascii:",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(float);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(float),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

off_t sf_bytes (sf_file file)
/*< Count the file data size (in bytes) >*/
{
    int st;
    off_t size;
    struct stat buf;
    
    if (0 == strcmp(file->dataname,"stdin")) return ((off_t) -1);

    if (NULL == file->dataname) {
	st = fstat(fileno(file->stream),&buf);
    } else {
	st = stat(file->dataname,&buf);
    }
    if (0 != st) sf_error ("%s: cannot find file size:",__FILE__);
    size = buf.st_size;
    return size;
}

off_t sf_tell (sf_file file)
/*< Find position in file >*/
{
    extern off_t ftello (FILE *stream);
    return ftello(file->stream);
}

FILE *sf_tempfile(char** dataname, const char* mode)
/*< Create a temporary file with a unique name >*/
{
    FILE *tmp;
    char *path;
    extern FILE * fdopen (int fd, const char *mode);
    extern int mkstemp (char *template);
    
    path = gettmpdatapath();
    if (NULL == path) path = getdatapath();
    *dataname = sf_charalloc (NAME_MAX+1);
    snprintf(*dataname,NAME_MAX,"%s%sXXXXXX",path,sf_getprog());
    tmp = fdopen(mkstemp(*dataname),mode);
    if (NULL == tmp) sf_error ("%s: cannot open %s:",__FILE__,*dataname);
    
    return tmp;
}

void sf_seek (sf_file file, off_t offset, int whence)
/*< Seek to a position in file. Follows fseek convention. >*/
{
    extern int fseeko(FILE *stream, off_t offset, int whence);

    if (0 > fseeko(file->stream,offset,whence))
	sf_error ("%s: seek problem:",__FILE__);
}

void sf_unpipe (sf_file file, off_t size) 
/*< Redirect a pipe input to a direct access file >*/
{
    off_t nbuf, len, bufsiz;
    char *dataname=NULL;
    FILE* tmp;
    char *buf;

    if (!(file->pipe)) return;
	
    tmp = sf_tempfile(&dataname,"wb");

    bufsiz = sf_bufsiz(file);
    buf = sf_charalloc(SF_MIN(bufsiz,size));

    while (size > 0) {
	nbuf = (bufsiz < size)? bufsiz : size;
	if (nbuf != fread(buf,1,nbuf,file->stream) ||
	    nbuf != fwrite(buf,1,nbuf,tmp))
	    sf_error ("%s: trouble unpiping:",__FILE__);
	size -= nbuf;
    }

    free(buf);

    if (NULL != file->dataname ) {
	len = strlen(dataname)+1;
	file->dataname = (char*) sf_realloc(file->dataname,len,sizeof(char));
	memcpy(file->dataname,dataname,len);
    }

    /*
      if (unlink(file->dataname))
      sf_warning ("%s: trouble removing %s:",__FILE__,file->dataname);
    */

    (void) fclose(file->stream);
    file->stream = freopen(dataname,"rb",tmp);

    if (NULL == file->stream)
	sf_error ("%s: Trouble reading data file %s:",__FILE__,dataname);
} 

void sf_close(void)
/*< Remove temporary file, created by sf_unpipe >*/
{
    if (NULL == infile || NULL == infile->dataname || !(infile->pipe)) return;
    
    if (0 != strcmp("stdin",infile->dataname) && 
	0 != unlink(infile->dataname))
	sf_warning ("%s: trouble removing %s:",__FILE__,infile->dataname);
}

/* 	$Id$	 */

