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
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sun) || defined (__sgi) || defined(hpux) || defined(__hpux)
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

typedef enum {SF_UCHAR, SF_CHAR, SF_INT, SF_FLOAT, SF_COMPLEX, SF_SHORT, SF_DOUBLE} sf_datatype;
typedef enum {SF_ASCII, SF_XDR, SF_NATIVE} sf_dataform;
/*^*/

#endif

struct sf_File {
    FILE *stream, *head; 
    char *dataname, *buf, *headname;
    sf_simtab pars;
    XDR xdr;
    enum xdr_op op;
    sf_datatype type;
    sf_dataform form;
    bool pipe, rw, dryrun;
};

/*@null@*/ static sf_file *infiles = NULL;
static size_t nfile=0, ifile=0;
static const int tabsize=10;
/*@null@*/ static char *aformat = NULL;
static size_t aline=8;
static bool error=true;

static bool getfilename (FILE *fd, char *filename);
static char* getdatapath (void);
static char* gettmpdatapath (void);
static bool readpathfile (const char* filename, char* datapath);
static void sf_input_error(sf_file file, const char* message, const char* name);

void sf_file_error(bool err)
/*< set error on opening files >*/
{
    error = err;
}

static void sf_input_error(sf_file file, const char* message, const char* name)
{
    if (error) 
	sf_error ("%s: %s %s:",__FILE__,message,name);
    sf_fileclose(file);
}

sf_file sf_input (/*@null@*/ const char* tag)
/*< Create an input file structure >*/
{
    int esize;
    sf_file file;
    char *filename, *format;
    size_t len;
    extern off_t ftello (FILE *stream);
	
    file = (sf_file) sf_alloc(1,sizeof(*file));
    file->dataname = NULL;
    
    if (NULL == tag || 0 == strcmp(tag,"in")) {
	file->stream = stdin;
	filename = NULL;
    } else {
	filename = sf_getstring (tag);
	if (NULL == filename) {
	    /* this option allows you to call function with 
	       sf_input("mydir/myinput.rsf");  Karl  */
	    len = strlen(tag)+1;
	    filename = sf_charalloc(len);
	    /* should change next line to strcpy or strncpy Karl */
	    memcpy(filename,tag,len);
	}
		
	file->stream = fopen(filename,"r");
	if (NULL == file->stream) {
	    sf_input_error(file,"Cannot read input (header) file",filename);
	    return NULL;
	}
    }
	
    file->buf = NULL;
    /*    setbuf(file->stream,file->buf); */
	
    /* create a parameter table */
    file->pars = sf_simtab_init (tabsize);
    file->head = sf_tempfile(&file->headname,"w+");
	
    /* read the parameter table from the file */
    sf_simtab_input (file->pars,file->stream,file->head);
	
    if (NULL == infiles) {
	infiles = (sf_file *) sf_alloc(1,sizeof(sf_file));
	infiles[0] = NULL;
	nfile=1;
    }
	
    if (NULL == filename) {
	infiles[0] = file;
    } else {
	free (filename);
	ifile++;
	if (ifile >= nfile) {
	    /* grow array */
	    nfile *= 2;
	    infiles = (sf_file *) realloc(infiles, nfile * sizeof(sf_file));
	    if (NULL == infiles) sf_error("%s: reallocation error",__FILE__);
	}
	infiles[ifile] = file;
    }
	
    filename = sf_histstring(file,"in");
    if (NULL == filename) {
    	sf_input_error (file,"No in= in file",tag);
	return NULL;
    }
    len = strlen(filename)+1;
    file->dataname = sf_charalloc(len);
    memcpy(file->dataname,filename,len);
	
    if (0 != strcmp(filename,"stdin")) {
	file->stream = freopen(filename,"rb",file->stream);
	if (NULL == file->stream) {
	    sf_input_error(file,"Cannot read data file",filename);
	    return NULL;
	}
    }
    free (filename);
	
    file->pipe = (bool) (-1 == ftello(file->stream));
    if (file->pipe && ESPIPE != errno) 
	sf_error ("%s: pipe problem:",__FILE__);
	
    file->op = XDR_DECODE;
	
    format = sf_histstring(file,"data_format");
    if (NULL == format) {
	if (!sf_histint(file,"esize",&esize) || 0 != esize) {
	    sf_input_error (file,"Unknown format in",tag);
	    return NULL;
	}
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
  Should do output after the first call to sf_input. >*/
{
    sf_file file;
    char *headname, *dataname, *path, *name, *format;
    size_t namelen;
    extern int mkstemp (char *tmpl);
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
        {
            free(file);
	    sf_error ("%s: Cannot write to header file %s:",__FILE__,headname);
        }
    }
	
    file->buf = NULL;
    /*    setbuf(file->stream,file->buf); */
	
    file->pars = sf_simtab_init (tabsize);
    file->head = NULL;
    file->headname = NULL;
	
    file->pipe = (bool) (-1 == ftello(file->stream));
    if (file->pipe && ESPIPE != errno) 
    {
        free(file);
	sf_error ("%s: pipe problem:",__FILE__);
    }
 
    dataname = sf_getstring("out");
    if (NULL == dataname)
	dataname = sf_getstring("--out");
	
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
	  if(0==strcmp(name,"/dev/null")){
	      file->dataname = sf_charalloc (7);
	      memcpy(file->dataname,"stdout",7);
	    } else {
	      namelen = strlen(file->dataname);
	      file->dataname[namelen]='@';
	      file->dataname[namelen+1]='\0';
	    }
	} else { /* Invent a name */
	  /* stdout is not a pipe, not /dev/null, not a file in this 
	     directory. 
	     One way to get here is to redirect io to a file not in this 
	     directory.  For example >../myfile.rsf.  In this case getfilename
	     cannot find the file from file->stream by looking in the current 
	     directory.  The function mkstemp is used to create a unique name 
	     to contain the binary data. */
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
	
    if (NULL == infiles) {
	infiles = (sf_file *) sf_alloc(1,sizeof(sf_file));
	infiles[0] = NULL;
	nfile=1;
    } 
	
    if (NULL != infiles[0]) { 
	if (NULL == infiles[0]->pars) sf_error("The input file was closed prematurely.");
	if (NULL != (format = sf_histstring(infiles[0],"data_format"))) {
	    sf_setformat(file,format);
	    free (format);
	}
    } else {
	sf_setformat(file,"native_float");
    }
	
    if (NULL != headname) free(headname);
	
    if (!sf_getbool("--readwrite",&(file->rw))) file->rw=false;
    if (!sf_getbool("--dryrun",&(file->dryrun))) file->dryrun=false;
	
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
        case SF_SHORT:
            return sizeof(short);
            break;
        case SF_DOUBLE:
            return sizeof(double);
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
    if (NULL != file->dataname) sf_putint(file,"esize",(int) sf_esize(file));
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

void sf_expandpars (sf_file file)
/*< add parameters from the command line >*/
{
    sf_simtab_expand(file->pars,sf_getpars());
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
	    if (NULL != file->dataname) sf_putint(file,"esize",0); /* for compatibility with SEPlib */
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
    } else if (NULL != strstr(format,"short")) {
	sf_settype(file,SF_SHORT);
    } else if (NULL != strstr(format,"double")) {
	sf_settype(file,SF_DOUBLE);
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
    struct stat buf_fstat_dev_null;
    FILE* fd_dev_null;
	
    dir = opendir(".");
    if (NULL == dir) sf_error ("%s: cannot open current directory:",__FILE__);
    
    if(0 > fstat(fileno(fp),&buf)) sf_error ("%s: cannot run fstat:",__FILE__);
    success = false;

    fd_dev_null=fopen("/dev/null","w");
    fstat(fileno(fd_dev_null),&buf_fstat_dev_null);
    if(buf_fstat_dev_null.st_dev==buf.st_dev){
	strcpy(filename,"/dev/null");
	success = true;
    } else {
      while (NULL != (dirp = readdir(dir))) {
	if (dirp->d_ino == buf.st_ino) { /* non-portable */
	  strcpy(filename,dirp->d_name);
	  success = true;
	  break;
	}
      }
    }
    fclose(fd_dev_null);
	
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
	
    free(path);
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

   13245 is more reasonable (Zhonghuan Chen)
   seems like option 5 would result in:
      if you move the header to another directory then you must also move 
      the data to the same directory. Confusing consequence. 
      (Karl Schleicher 5/28/2013) 
*/
{
    char *path, *penv, *home, file[PATH_MAX];
    
    path = sf_getstring ("datapath");
    if (NULL != path) return path;
	
    path = sf_charalloc(PATH_MAX+1);
	
    if (readpathfile (".datapath",path)) return path;
    
    penv = getenv("DATAPATH");
    if (NULL != penv) {
	strncpy(path,penv,PATH_MAX);
	return path;
    }
	
    home = getenv("HOME");
    if (NULL != home) {
	(void) snprintf(file,PATH_MAX,"%s/.datapath",home);
	if (readpathfile (file,path)) return path;
    }
	
    /* could next 4 line be simplified to:
       return "./"; ?  or just return NULL;??     Karl Schleicher 5/28/2013 */
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
    if (NULL == file) return;
    
    if (file->stream != stdin && 
	file->stream != stdout && 
	file->stream != NULL) {
	(void) fclose (file->stream);
	file->stream = NULL;
    }
	
    if (file->head != NULL) {
	(void) unlink (file->headname);
	(void) fclose (file->head);
	file->head = NULL;
	free(file->headname);
	file->headname = NULL;
    }
	
    if (NULL != file->pars) {
	sf_simtab_close (file->pars);
	file->pars = NULL;
    }
	
    if (NULL != file->buf) {
	free (file->buf);
	file->buf = NULL;
    }
	
    if (NULL != file->dataname) {
	free (file->dataname);
	file->dataname = NULL;
    }
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

bool sf_histlargeint (sf_file file, const char* key,/*@out@*/ off_t* par)
/*< read a sf_largeint parameter from file >*/ 
{
    return sf_simtab_getlargeint (file->pars,key,par);
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
    if (0 >= fprintf(file->stream,"%s\t%s\t%s:\t%s@%s\t%s\n",
		     RSF_VERSION,
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
        case SF_SHORT:
	    switch (file->form) {
		case SF_ASCII:
		    sf_putstring(file,"data_format","ascii_short");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_short");
		    break;
		default:
		    sf_putstring(file,"data_format","native_short");
		    break;
	    }
	    break;
        case SF_DOUBLE:
	    switch (file->form) {
		case SF_ASCII:
		    sf_putstring(file,"data_format","ascii_double");
		    break;
		case SF_XDR:
		    sf_putstring(file,"data_format","xdr_double");
		    break;
		default:
		    sf_putstring(file,"data_format","native_double");
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
    (void) fflush(file->stream);
	
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
	
    if (NULL == file->dataname) 
	sf_warning("%s: putint to a closed file",__FILE__);
    snprintf(val,256,"%d",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putints (sf_file file, const char* key, const int* par, size_t n)
/*< put an int array of size n parameter to a file >*/
{
    int i;
    char val[1024], *v;
	
    if (NULL == file->dataname) 
	sf_warning("%s: putints to a closed file",__FILE__);
    v = val;
    for (i=0; i < (int) n-1; i++) {
	v += snprintf(v,1024,"%d,",par[i]);
    }
    snprintf(v,1024,"%d",par[n-1]);
	
    sf_simtab_enter (file->pars,key,val);
}

void sf_putlargeint (sf_file file, const char* key, off_t par)
/*< put a sf_largeint parameter to a file >*/
{
    char val[256];
	
    if (NULL == file->dataname) 
	sf_warning("%s: putflargeint to a closed file",__FILE__);
    snprintf(val,256,"%ld",(long int) par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putfloat (sf_file file, const char* key,float par)
/*< put a float parameter to a file >*/
{
    char val[256];
	
    if (NULL == file->dataname) 
	sf_warning("%s: putfloat to a closed file",__FILE__);
    snprintf(val,256,"%g",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putstring (sf_file file, const char* key,const char* par)
/*< put a string parameter to a file >*/
{
    char *val;
    
    if (NULL == file->dataname) 
	sf_warning("%s: putstring to a closed file",__FILE__);
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
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
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
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
	    break;
    }
}

void sf_charwrite (char* arr, size_t size, sf_file file)
/*< write a char array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
	
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
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
	
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
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
	    break;
    }
}

int sf_try_charread(const char* test, sf_file file)
/*< check if you can read test word >*/
{
    int size, got, cmp;
    char* arr;

    size = strlen(test);
    arr = sf_charalloc(size+1);
    got = fread(arr,sizeof(char),size+1,file->stream);
    if (got != size) return 1;
    arr[size] = '\0';
    cmp = strncmp(arr,test,size);
    free(arr);

    return cmp;
}

int sf_try_charread2 (/*@out@*/ char* arr, size_t size, sf_file file)
/*< try to read size bytes.  return number bytes read >*/
{
    return fread(arr,sizeof(char),size,file->stream);
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
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
	    break;
    }
}

void sf_intwrite (int* arr, size_t size, sf_file file)
/*< write an int array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
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
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
	    break;
    }
}

void sf_shortread (/*@out@*/ short* arr, size_t size, sf_file file)
/*< read a short array arr[size] from file >*/
{
    char* buf;
    size_t i, left, nbuf, got, bufsiz;
	
    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		if (EOF==fscanf(file->stream,"%hd",arr+i))
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
				nbuf/sizeof(short),sizeof(short),
				(xdrproc_t) xdr_int))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(short),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
	    break;
    }
}

void sf_shortwrite (short* arr, size_t size, sf_file file)
/*< write a short array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
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
	    size *= sizeof(short);
	    buf = (char*)arr+size;
	    bufsiz = sf_bufsiz(file);
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (bufsiz < left)? bufsiz : left;
		(void) xdr_setpos(&(file->xdr),0);
		if (!xdr_vector(&(file->xdr),buf-left,
				nbuf/sizeof(short),sizeof(short),
				(xdrproc_t) xdr_int))
		    sf_error ("sf_file: trouble writing xdr");
		fwrite(file->buf,1,nbuf,file->stream);
	    }
	    break;
	default:
	    fwrite(arr,sizeof(short),size,file->stream);
	    break;
    }
}

void sf_floatwrite (float* arr, size_t size, sf_file file)
/*< write a float array arr[size] to file >*/
{
    char* buf;
    size_t i, left, nbuf, bufsiz;
	
    if (NULL != file->dataname) sf_fileflush (file,infiles[0]);
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
		sf_error ("%s: trouble reading: %lu of %lu",__FILE__,got,size);
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
    int stemp;
    extern FILE * fdopen (int fd, const char *mode);
    extern int mkstemp (char *tmpl);
    
    path = gettmpdatapath();
    if (NULL == path) path = getdatapath();
    *dataname = sf_charalloc (NAME_MAX+1);
    snprintf(*dataname,NAME_MAX,"%s%sXXXXXX",path,sf_getprog());

    free(path);
    
    stemp = mkstemp(*dataname);
    if (stemp < 0) sf_error ("%s: cannot create %s:",__FILE__,*dataname);
	
    tmp = fdopen(stemp,mode);
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

FILE* sf_filestream (sf_file file)
/*< Returns file descriptor to a stream >*/
{
    return file ? file->stream : NULL;
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
/*< Remove temporary files >*/
{
    int i;
    sf_file file;
	
    if (NULL == infiles) return;
    
    for (i=0; i <= (int) ifile; i++) {
	file = infiles[i];
		
	if (NULL != file && 
	    NULL != file->dataname && 
	    file->pipe &&
	    0 != strcmp("stdin",file->dataname) && 
	    0 != unlink(file->dataname))
	    sf_warning ("%s: trouble removing %s:",__FILE__,file->dataname);
	    
	sf_fileclose(file);
	free(file);
    }
    free(infiles);
    infiles=NULL;
    ifile=nfile=0;
}

/* 	$Id$	 */

