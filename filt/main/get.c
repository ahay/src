/* Output parameters from the header.

Takes: < file.rsf par1 par2 ...
*/

#include <stdio.h>
#include <string.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int i;
    char *string, *key;
    bool parform;
    sf_file in;

    sf_init (argc,argv);
    in = sf_input ("in");
    
    if(!sf_getbool("parform",&parform)) parform=true;
    /* If y, print out parameter=value. If n, print out value. */

    for (i = 1; i < argc; i++) {
	key = argv[i];
	if (NULL != strchr(key,'=')) continue;
	string = sf_histstring(in,key);
	if (NULL == string) {
	    sf_warning("No key %s",key);
	} else {
	    if (parform) printf ("%s=",key);
	    printf("%s\n",string);
	} 
    }
    
    exit(0);
}

/* 	$Id: get.c,v 1.2 2003/09/29 14:34:56 fomels Exp $	 */

