#include <stdlib.h>
#include <string.h>

#include <sys/utsname.h>
#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>

#include "getpar.h"
#include "c99.h"
#include "simtab.h"
#include "alloc.h"
#include "error.h"

static sf_simtab pars;
static char *prog = NULL;
static char *user = NULL;
static char *host = NULL;

void sf_init(int argc,char *argv[]) {
    int ic;
    struct passwd* pass;
    struct utsname uhost;
    FILE *fp;
    size_t len;

    pars = sf_simtab_init (argc);

    /* set prog */
    prog = strrchr(argv[0],'/');
    prog = (NULL == prog)? argv[0]:prog+1;

    /* set user */
    user = getlogin();
    if (NULL == user) {
	pass = getpwuid(geteuid());
	if (NULL == pass) {
	    user = sf_charalloc(1);
	    user[0] = '\0';
	} else {
	    user = pass->pw_name;
	}
    }

    /* set host */
    if (0 > uname (&uhost))
	sf_error ("%s: Cannot get hostname:",__FILE__);
    len = strlen(uhost.nodename)+1;
    host = sf_charalloc(len);
    memcpy(host,uhost.nodename,len);

    for (ic=1; ic < argc; ic++) {
	if (0 == strncmp(argv[ic],"par=",4)) {
	    fp = fopen(argv[ic]+4,"r");
	    if (NULL == fp) 
		sf_error ("%s: Cannot open par file %s:",__FILE__,argv[ic]+4);
	    sf_simtab_input(pars,fp);
	    (void) fclose (fp);
	} else {
	    sf_simtab_put(pars,argv[ic]);
	}
    }
}

void sf_close (void)
{
    sf_simtab_close (pars);
}

char* sf_getprog (void) {
    return prog;
}

char* sf_getuser (void) {
    return user;
}

char* sf_gethost (void) {
    return host;
}

bool sf_getint (const char* key,/*@out@*/ int* par) {
    return sf_simtab_getint (pars,key,par);
}

bool sf_getints (const char* key,/*@out@*/ int* par,size_t n) {
    return sf_simtab_getints(pars,key,par,n);
} 

bool sf_getfloat (const char* key,/*@out@*/ float* par) {
    return sf_simtab_getfloat (pars,key,par);
}

bool sf_getfloats (const char* key,/*@out@*/ float* par,size_t n) {
    return sf_simtab_getfloats (pars,key,par,n);
}

char* sf_getstring (const char* key) {
    return sf_simtab_getstring(pars,key);
}

bool sf_getstrings (const char* key,/*@out@*/ char** par,size_t n) {
    return sf_simtab_getstrings(pars,key,par,n);
} 

bool sf_getbool (const char* key,/*@out@*/ bool* par) {
    return sf_simtab_getbool(pars,key,par);
} 
    
bool sf_getbools (const char* key,/*@out@*/ bool* par,size_t n) {
    return sf_simtab_getbools(pars,key,par,n);
} 
 


