#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "files.h"
#include "c99.h"
#include "file.h"
#include "simtab.h"
#include "error.h"
#include "alloc.h"

int sf_filedims (sf_file file, /*@out@*/ int *n) 
{
    int i, dim;
    char key[3];

    dim = 1;
    for (i=0; i < SF_MAX_DIM; i++) {
	snprintf(key,3,"n%d",i+1);
	if (!sf_histint(file,key,n+i)) break;
	if (n[i] > 1) dim=i+1;
    }
    return dim;
}

int sf_filesize (sf_file file) 
{    
    return sf_leftsize (file, 0);
}

int sf_leftsize (sf_file file, int dim) 
{
    int size, ni;
    char key[3];

    for (size=1; dim < SF_MAX_DIM; dim++, size *= ni) {
	snprintf(key,3,"n%d",dim+1);
	if (!sf_histint(file,key,&ni)) break;
    }
    return size;
}

void sf_cp(sf_file in, sf_file out)
{
    int esize;
    long nsiz, nbuf;
    char buf[BUFSIZ];
    
    nsiz = sf_bytes (in);
    if (nsiz < 0) { /* reading from "stdin" */
	nsiz = sf_filesize (in);
	if (!sf_histint(in,"esize",&esize) || esize <= 0)
	    sf_error("%s: need esize > 0 in input",__FILE__);
	nsiz *= esize;
    }

    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    for (nbuf = BUFSIZ; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;
	sf_read (buf,1,nbuf,in);
	sf_write (buf,1,nbuf,out);
    }
}

void sf_rm(const char* filename, bool force, bool verb, bool inquire)
{
    int c, c2;
    char cc, *in;
    FILE *file, *query;
    sf_simtab tab;
    struct stat buf;
    mode_t mod;
    const int tabsize=10;
    
    tab = sf_simtab_init (tabsize);
    query = fopen ("/dev/tty","w+");
    if (inquire) {
	if (NULL == query) sf_error ("%s: Cannot open terminal",__FILE__);
	setbuf (query,NULL);
	fprintf (query,"sf_rm: Remove '%s'? ",filename);
	c2 = c = getc (query);
	while (c2 != EOF && (cc = (char) c2) != '\n' && cc != '\0') 
	    c2 = getc(query);
	cc = (char) c;
	if ('y' != cc && 'Y' != cc) return;
    }
    if (verb) sf_warning("sf_rm: Removing header %s",filename);
    file = fopen (filename,"r");
    if (NULL == file) sf_error ("%s: Cannot open file %s:",__FILE__,filename);
    sf_simtab_input (tab,file);
    (void) fclose (file);
    in = sf_simtab_getstring (tab,"in");
    if (NULL == in) sf_error ("%s:  File %s has no in=",__FILE__,filename);
    if (0 != remove(filename)) 
	sf_error ("%s: Trouble removing header file %s:",__FILE__,filename);
	    
    if (verb) sf_warning("sf_rm: Removing data %s",in);
    if (!force) {
	if (0 != stat(in,&buf)) 
	    sf_error ("%s: Trouble with file %s:",__FILE__,in);
	/* (owner and can write) or (not-owner and others can write) */
	mod = (buf.st_uid == getuid())? S_IWUSR:S_IWOTH;
	if (0 == (buf.st_mode & mod)) {		    
	    fprintf (query,"sf_rm: Remove protected file '%s'? ",in);
	    c2 = c = getc (query);
	    while (c2 != EOF && (cc = (char) c2) != '\n' && cc != '\0') 
		c2 = getc(query);
	    cc = (char) c;
	    if ('y' != cc && 'Y' != cc) return;
	}
    }
    if (0 != remove(in)) 
	sf_error ("%s: Trouble removing data file %s:",__FILE__,in);
    sf_simtab_close (tab);
}

/* 	$Id: files.c,v 1.2 2003/09/29 14:34:55 fomels Exp $	 */
