#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __GNUC__
#ifndef alloca
#define alloca __builtin_alloca
#endif
#else /* not GNU C  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#endif
#endif

#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>

#include <rpc/types.h>
/* #include <rpc/rpc.h> */
#include <rpc/xdr.h>

#include "file.h"
#include "c99.h"
#include "getpar.h"
#include "simtab.h"
#include "alloc.h"
#include "error.h"

/* BSD - MAXNAMELEN, Posix - NAME_MAX */
#ifndef NAME_MAX
#ifdef MAXNAMELEN
#define	NAME_MAX MAXNAMELEN
#endif
#ifdef FILENAME_MAX
#define NAME_MAX FILENAME_MAX
#endif
#endif

struct sf_File {
    FILE *stream; 
    char *dataname, *buf;
    sf_simtab pars;
    XDR *xdr;
    enum xdr_op op;
    sf_datatype type;
    sf_dataform form;
    bool pipe;
};

/*@null@*/ static sf_file infile=NULL;
static const int tabsize=10;
/*@null@*/ static char *aformat = NULL;
static size_t aline=8;

static bool getfilename (FILE *fd, char *filename);
static char* getdatapath (void);
static bool readpathfile (const char* filename, char* datapath);

sf_file sf_input (/*@null@*/ const char* tag)
{
    int esize;
    sf_file file;
    char *filename, *format;
    size_t len;

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

    file->pars = sf_simtab_init (tabsize);
    sf_simtab_input (file->pars,file->stream);
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

    file->pipe = (-1L == fseek(file->stream,0L,SEEK_CUR));
    if (file->pipe && ESPIPE != errno) 
	sf_error ("%s: pipe problem:",__FILE__);

    file->op = XDR_DECODE;
    file->xdr = NULL;

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

/* Should do output after input */
sf_file sf_output (/*@null@*/ const char* tag)
{
    sf_file file;
    char *headname, *dataname, *path, *name, *format;
    size_t namelen;

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

    file->pipe = (-1L == fseek(file->stream,0L,SEEK_CUR));
    if (file->pipe && ESPIPE != errno) 
	sf_error ("%s: pipe problem:",__FILE__);

    dataname = sf_getstring("out");

    if (file->pipe) {
	file->dataname = sf_charalloc (7);
	memcpy(file->dataname,"stdout",7);
    } else if (NULL == dataname) {
	path = getdatapath();
	if (NULL == path) 
	    sf_error("%s: Cannot find datapath",__FILE__);
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
	    (void) unlink(file->dataname);
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
    file->xdr = NULL;

    if (NULL != infile && 
	NULL != (format = sf_histstring(infile,"data_format"))) {
	sf_setformat(file,format);
	free (format);
    } else {
	sf_setformat(file,"native_float");
    }

    if (NULL != headname) free(headname);
    return file;
}

sf_datatype sf_gettype (sf_file file)
{
    return file->type;
}

sf_dataform sf_getform (sf_file file)
{
    return file->form;
}

void sf_settype (sf_file file, sf_datatype type)
{
    file->type = type;
}

void sf_setform (sf_file file, sf_dataform form)
{
    file->form = form;
}

void sf_setformat (sf_file file, const char* format)
{
    if (NULL != strstr(format,"float")) {
	file->type = SF_FLOAT;
	sf_putint(file,"esize",(int) sizeof(float));
    } else if (NULL != strstr(format,"int")) {
	file->type = SF_INT;
	sf_putint(file,"esize",(int) sizeof(int));
    } else if (NULL != strstr(format,"complex")) {
	file->type = SF_COMPLEX;
	sf_putint(file,"esize",(int) sizeof(float complex));
    } else {
	file->type = SF_CHAR;
	sf_putint(file,"esize",(int) sizeof(char));
    }
    
    if (0 == strncmp(format,"ascii_",6)) {
	file->form = SF_ASCII;
	if (NULL != file->xdr) {
	    free (file->xdr);
	    free (file->buf);
	}
	file->xdr = NULL;
	sf_putint(file,"esize",0);
    } else if (0 == strncmp(format,"xdr_",4)) {
	file->form = SF_XDR;
	if (NULL == file->xdr) {
	    file->xdr = (XDR*) sf_alloc(1,sizeof(XDR));
	    file->buf = sf_charalloc(BUFSIZ);
	    xdrmem_create(file->xdr,file->buf,BUFSIZ,file->op);
	}
    } else {
	file->form = SF_NATIVE;
	if (NULL != file->xdr) {
	    free (file->xdr);
	    free (file->buf);
	}
	file->xdr = NULL;
    }
}
    
static bool getfilename (FILE* fp, char *filename)
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

    return success;
}

/* datapath rules:
   1. check datapath= on the command line
   2. check DATAPATH environmental variable
   3. check .datapath file in the current directory
   4. check .datapath in the home directory
*/

static char* getdatapath (void) 
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
	snprintf(file,PATH_MAX,"%s/.datapath",home);
	if (readpathfile (file,path)) return path;
    }

    return NULL;
}

static bool readpathfile (const char* filename, char* datapath) 
{
    FILE *fp;
    char format[PATH_MAX];

    fp = fopen(filename,"r");
    if (NULL == fp) return false;

    if (0 >= fscanf(fp,"datapath=%s",datapath))
	sf_error ("No datapath found in file %s",filename);

    snprintf(format,PATH_MAX,"%s datapath=%%s",sf_gethost());
    (void) fscanf(fp,format,datapath);

    (void) fclose (fp);
    return true;
}

void sf_fileclose (sf_file file) 
{
    if (file->stream != stdin && 
	file->stream != stdout && 
	file->stream != NULL) (void) fclose (file->stream);
    if (NULL != file->pars) sf_simtab_close (file->pars);
    if (NULL != file->xdr) {
	free (file->xdr);
	free (file->buf);
    }
    if (NULL != file->dataname) free (file->dataname);
    free (file);
}

bool sf_histint (sf_file file, const char* key,/*@out@*/ int* par) 
{
    return sf_simtab_getint (file->pars,key,par);
}

bool sf_histints (sf_file file, const char* key,/*@out@*/ int* par,size_t n) 
{
    return sf_simtab_getints (file->pars,key,par, n);
}

bool sf_histfloat (sf_file file, const char* key,/*@out@*/ float* par) 
{
    return sf_simtab_getfloat (file->pars,key,par);
}

bool sf_histfloats (sf_file file, const char* key,
		    /*@out@*/ float* par,size_t n) 
{
    return sf_simtab_getfloats (file->pars,key,par, n);
}

bool sf_histbool (sf_file file, const char* key,/*@out@*/ bool* par) 
{
    return sf_simtab_getbool (file->pars,key,par);
}

bool sf_histbools (sf_file file, const char* key,
		   /*@out@*/ bool* par, size_t n) 
{
    return sf_simtab_getbools (file->pars,key,par, n);
}

char* sf_histstring (sf_file file, const char* key) 
{
    return sf_simtab_getstring (file->pars,key);
}

void sf_fileflush (sf_file file, sf_file src)
{
#if 0
    int i, n;
    float f;
    char key[8], *val;
#endif
    time_t tm;
 
    if (NULL == file->dataname) return;
    
    tm = time (NULL);
    if (0 >= fprintf(file->stream,"%s:\t%s@%s\t%s\n",
		     sf_getprog(),
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
	default:
	    sf_putstring(file,"data_format",
			 (NULL==file->xdr)? "native_byte":"xdr_byte");
	    break;
    }    
    if (NULL != src && NULL != src->pars)
	sf_simtab_output(src->pars,file->stream);
#if 0
        for (i=1; i <= SF_MAX_DIM; i++) {
	    snprintf(key,4,"n%d",i);
	    if (!sf_simtab_getint(src->pars,key,&n)) break;
	    if (0 >= fprintf(file->stream,"\t%s=%d\n",key,n))
		sf_error ("%s: cannot flush %s:",__FILE__,key);
	    snprintf(key,4,"o%d",i);
	    if (sf_simtab_getfloat(src->pars,key,&f) &&
		0 >= fprintf(file->stream,"\t%s=%g\n",key,f))
		sf_error ("%s: cannot flush %s:",__FILE__,key);
	    snprintf(key,4,"d%d",i);
	    if (sf_simtab_getfloat(src->pars,key,&f) &&
		0 >= fprintf(file->stream,"\t%s=%g\n",key,f))
		sf_error ("%s: cannot flush %s:",__FILE__,key);
	    snprintf(key,8,"label%d",i);
	    if (NULL != (val=sf_simtab_getstring(src->pars,key)) &&
		0 >= fprintf(file->stream,"\t%s=\"%s\"\n",key,val))
		sf_error ("%s: cannot flush %s:",__FILE__,key);
	}
#endif
    
    sf_simtab_output(file->pars,file->stream);
   
    if (0==strcmp(file->dataname,"stdout")) { 
        /* keep stream, write the header end code */
	fprintf(file->stream,"\tin=\"stdin\"\n\n%c%c%c",
		SF_EOL,SF_EOL,SF_EOT);
    } else {
	file->stream = freopen(file->dataname,"wb",file->stream);
	if (NULL == file->stream) 
	    sf_error ("%s: Cannot write to data file %s:",
		      __FILE__,file->dataname);	
    }

    free (file->dataname);
    file->dataname=NULL;
}

void sf_putint (sf_file file, const char* key, int par)
{
    char val[256];

    snprintf(val,256,"%d",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putints (sf_file file, const char* key, const int* par, size_t n)
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


void sf_putfloat (sf_file file, const char* key,float par)
{
    char val[256];

    snprintf(val,256,"%g",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putstring (sf_file file, const char* key,const char* par)
{
    char *val;
    
    val = (char*) alloca(strlen(par)+3); 
    sprintf(val,"\"%s\"",par);
    sf_simtab_enter (file->pars,key,val);
}

void sf_putline (sf_file file, const char* line)
{
    if (0 >= fprintf(file->stream,"\t%s\n",line))
	sf_error ("%s: cannot put line '%s':",__FILE__,line);
}

void sf_setaformat (const char* format, int line)
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

/*
void sf_write (void* arr, size_t esize, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf;
    bool_t success;
    float complex c;

    if (NULL != file->dataname) sf_fileflush (file,infile);
    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		switch(file->type) {
		    case SF_INT:
			for (i=size-left; i < size-left+nbuf; i++) {
			    if (EOF==fprintf(file->stream,
					     (NULL != aformat)? aformat:"%d ",
					     ((int*)arr)[i]))
				sf_error ("%s: trouble writing ascii",__FILE__);
			}
			break;
		    case SF_FLOAT:
			for (i=size-left; i < size-left+nbuf; i++) {
			    if (EOF==fprintf(file->stream,
					     (NULL != aformat)? aformat:"%g ",
					     ((float*)arr)[i]))
				sf_error ("%s: trouble writing ascii",__FILE__);
			}
			break;
		    case SF_COMPLEX:
			for (i=size-left; i < size-left+nbuf; i++) {
			    c = ((float complex*)arr)[i];
			    if (EOF==fprintf(file->stream,
					     (NULL != aformat)? 
					     aformat:"(%g,%g) ",
					     crealf(c),cimagf(c)))
				sf_error ("%s: trouble writing ascii",__FILE__);
			}
			break;
		    default:
			for (i=size-left; i < size-left+nbuf; i++) {
			    if (EOF==fputc(((char*)arr)[i],file->stream))
				sf_error ("%s: trouble writing ascii",__FILE__);
			}
			break;
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= esize;
	    buf = (char*)arr+size;
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		switch(file->type) {
		    case SF_INT:
			success=xdr_vector(file->xdr,buf-left,
				     nbuf/sizeof(int),sizeof(int),
				     (xdrproc_t) xdr_int);
			break;
		    case SF_FLOAT:
		    case SF_COMPLEX:
			success=xdr_vector(file->xdr,buf-left,
				nbuf/sizeof(float),sizeof(float),
					(xdrproc_t) xdr_float);
			break;
		    default:
			success=xdr_opaque(file->xdr,buf-left,nbuf);
			break;
		}       
		if (0 == success) sf_error ("sf_file: trouble writing xdr");
		if (nbuf != fwrite(file->buf,1,nbuf,file->stream)) 
		    sf_error ("%s: trouble writing:",__FILE__);
	    }
	    break;
	default:
	    if (size != fwrite(arr,esize,size,file->stream)) 
		sf_error ("%s: trouble writing:",__FILE__);
	    break;
    }
}

void sf_read (void* arr, size_t esize, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf, got;
    bool_t success;
    int c;
    float re, im;

    switch (file->form) {
	case SF_ASCII:
	    switch (file->type) {
		case SF_INT:
		    for (i = 0; i < size; i++) {
			if (EOF==fscanf(file->stream,"%d",(int*)arr+i))
			    sf_error ("%s: trouble reading ascii:",__FILE__);
		    }
		    break;
		case SF_FLOAT:
		    for (i = 0; i < size; i++) {
			if (EOF==fscanf(file->stream,"%g",(float*)arr+i))
			    sf_error ("%s: trouble reading ascii:",__FILE__);
		    }
		    break;
		case SF_COMPLEX:
		    for (i = 0; i < size; i++) {
			if (EOF==fscanf(file->stream,"(%g,%g)",&re,&im))
			    sf_error ("%s: trouble reading ascii:",__FILE__);
			((float complex*)arr)[i]=re+I*im;
		    }
		    break;
		default:
		    for (i = 0; i < size; i++) {
			c=fgetc(file->stream);
			if (EOF==c)
			    sf_error ("%s: trouble reading ascii:",__FILE__);
			((char*)arr)[i]= (char) c;
		    }
		    break;
	    }
	    break;
	case SF_XDR:
	    size *= esize;
	    buf = (char*)arr+size;
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		switch (file->type) {
		    case SF_INT:
			success=xdr_vector(file->xdr,buf-left,
				     nbuf/sizeof(int),sizeof(int),
				     (xdrproc_t) xdr_int);
			break;
		    case SF_FLOAT:
		    case SF_COMPLEX:
			success=xdr_vector(file->xdr,buf-left,
				     nbuf/sizeof(float),sizeof(float),
				     (xdrproc_t) xdr_float);
			break;
		    default:
			success=xdr_opaque(file->xdr,buf-left,nbuf);
			break;
		}
		if (0==success) sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,esize,size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}
*/

void sf_complexwrite (float complex* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf;
    float complex c;

    if (NULL != file->dataname) sf_fileflush (file,infile);
    switch(file->form) {
	case SF_ASCII:
	    for (left = size; left > 0; left-=nbuf) {
		nbuf = (aline < left)? aline: left;
		for (i=size-left; i < size-left+nbuf; i++) {
		    c = arr[i];
		    if (EOF==fprintf(file->stream,
				     (NULL != aformat)? 
				     aformat:"(%g,%g) ",
				     crealf(c),cimagf(c)))
			sf_error ("%s: trouble writing ascii",__FILE__);
		}
		if (EOF==fprintf(file->stream,"\n"))
		    sf_error ("%s: trouble writing ascii",__FILE__);
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(float complex);
	    buf = (char*)arr+size;
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (!xdr_vector(file->xdr,buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("sf_file: trouble writing xdr");
		if (nbuf != fwrite(file->buf,1,nbuf,file->stream)) 
		    sf_error ("%s: trouble writing:",__FILE__);
	    }
	    break;
	default:
	    if (size != fwrite(arr,sizeof(float complex),size,file->stream)) 
		sf_error ("%s: trouble writing:",__FILE__);
	    break;
    }
}

void sf_complexread (/*@out@*/ float complex* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf, got;
    float re, im;

    switch (file->form) {
	case SF_ASCII:
	    for (i = 0; i < size; i++) {
		if (EOF==fscanf(file->stream,"(%g,%g)",&re,&im))
		    sf_error ("%s: trouble reading ascii:",__FILE__);
		arr[i]=re+I*im;
	    }
	    break;
	case SF_XDR:
	    size *= sizeof(float complex);
	    buf = (char*)arr+size;
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(file->xdr,buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("%s: trouble reading xdr",__FILE__);
	    }
	    break;
	default:
	    got = fread(arr,sizeof(float complex),size,file->stream);
	    if (got != size) 
		sf_error ("%s: trouble reading: %d of %d",__FILE__,got,size);
	    break;
    }
}

void sf_charwrite (char* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf;

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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (!xdr_opaque(file->xdr,buf-left,nbuf))
		    sf_error ("sf_file: trouble writing xdr");
		if (nbuf != fwrite(file->buf,1,nbuf,file->stream)) 
		    sf_error ("%s: trouble writing:",__FILE__);
	    }
	    break;
	default:
	    if (size != fwrite(arr,sizeof(char),size,file->stream)) 
		sf_error ("%s: trouble writing:",__FILE__);
	    break;
    }
}

void sf_charread (/*@out@*/ char* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf, got;
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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_opaque(file->xdr,buf-left,nbuf))
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

void sf_intwrite (int* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf;

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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (!xdr_vector(file->xdr,buf-left,
				nbuf/sizeof(int),sizeof(int),
				(xdrproc_t) xdr_int))
		    sf_error ("sf_file: trouble writing xdr");
		if (nbuf != fwrite(file->buf,1,nbuf,file->stream)) 
		    sf_error ("%s: trouble writing:",__FILE__);
	    }
	    break;
	default:
	    if (size != fwrite(arr,sizeof(int),size,file->stream)) 
		sf_error ("%s: trouble writing:",__FILE__);
	    break;
    }
}

void sf_intread (/*@out@*/ int* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf, got;

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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(file->xdr,buf-left,
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
{
    char* buf;
    size_t i, left, nbuf;

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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (!xdr_vector(file->xdr,buf-left,
				nbuf/sizeof(float),sizeof(float),
				(xdrproc_t) xdr_float))
		    sf_error ("sf_file: trouble writing xdr");
		if (nbuf != fwrite(file->buf,1,nbuf,file->stream)) 
		    sf_error ("%s: trouble writing:",__FILE__);
	    }
	    break;
	default:
	    if (size != fwrite(arr,sizeof(float),size,file->stream)) 
		sf_error ("%s: trouble writing:",__FILE__);
	    break;
    }
}

void sf_floatread (/*@out@*/ float* arr, size_t size, sf_file file)
{
    char* buf;
    size_t i, left, nbuf, got;

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
	    for (left = size; left > 0; left -= nbuf) {
		nbuf = (BUFSIZ < left)? BUFSIZ : left;
		(void) xdr_setpos(file->xdr,0);
		if (nbuf != fread(file->buf,1,nbuf,file->stream))
		    sf_error ("%s: trouble reading:",__FILE__);
		if (!xdr_vector(file->xdr,buf-left,
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

long sf_bytes (sf_file file)
{
    int st;
    long size;
    struct stat buf;
    
    if (0 == strcmp(file->dataname,"stdin")) return -1L;

    if (NULL == file->dataname) {
	st = fstat(fileno(file->stream),&buf);
    } else {
	st = stat(file->dataname,&buf);
    }
    if (0 != st) sf_error ("%s: cannot find file size:",__FILE__);
    size = buf.st_size;
    return size;
}
    
long sf_tell (sf_file file)
{
    return ftell(file->stream);
}

FILE *sf_tempfile(char** dataname)
{
    FILE *tmp;
    char *path;
    
    path = getdatapath();
    if (NULL == path) 
	sf_error ("%s: Cannot find datapath",__FILE__);
    *dataname = sf_charalloc (NAME_MAX+1);
    snprintf(*dataname,NAME_MAX,"%s%sXXXXXX",path,sf_getprog());
    tmp = fdopen(mkstemp(*dataname),"wb");
    if (NULL == tmp) sf_error ("%s: cannot open %s:",__FILE__,*dataname);

    return tmp;
}

void sf_seek (sf_file file, long offset, int whence)
{
    if (0 > fseek(file->stream,offset,whence))
	sf_error ("%s: seek problem:",__FILE__);
}

void sf_unpipe (sf_file file, size_t size) 
{
    size_t nbuf, len;
    char *dataname=NULL;
    FILE* tmp;
    char buf[BUFSIZ];

    if (!(file->pipe)) return;
	
    tmp = sf_tempfile(&dataname);

    while (size > 0) {
	nbuf = (BUFSIZ < size)? BUFSIZ : size;
	if (nbuf != fread(buf,1,nbuf,file->stream) ||
	    nbuf != fwrite(buf,1,nbuf,tmp))
	    sf_error ("%s: trouble unpiping:",__FILE__);
	size -= nbuf;
    }

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
{
    if (NULL == infile || NULL == infile->dataname || !(infile->pipe)) return;
    
    if (strcmp("stdin",infile->dataname) && unlink(infile->dataname))
	sf_warning ("%s: trouble removing %s:",__FILE__,infile->dataname);
}

/*
FILE* sf_direct (const sf_file file)
{
    return file->pipe? sf_tempfile(&(file->dataname)): file->stream;
}

void sf_pipe (sf_file file, FILE* tmp, size_t size) 
{
    char buf[BUFSIZ];
    size_t nbuf;

    if (tmp == file->stream) return;

    tmp = freopen(file->dataname,"rb",tmp);
    if (NULL == tmp)
	sf_error ("%s: Trouble reading data file %s:",__FILE__,file->dataname);

    while (size > 0) {
	nbuf = (BUFSIZ < size)? BUFSIZ : size;
	if (nbuf != fread(buf,1,nbuf,tmp) ||
	    nbuf != fwrite(buf,1,nbuf,file->stream))
	    sf_error ("%s: trouble piping from %s:",__FILE__,file->dataname);
	size -= nbuf;
    }

    (void) fclose(tmp);
}
*/

/* 	$Id$	 */

