/* Remove RSF files together with their data.

Takes: file1.rsf [file2.rsf ...] [-i] [-v] [-f] 

Mimics the standard Unix rm command.

See also: sfmv, sfcp.
*/

#include <string.h>

#include <rsf.h>

int main (int argc, char *argv[])
{
    int i;
    char *arg;
    bool force = false, verb = false, inquire = false;
   
    sf_init(argc,argv);

    for (i=1; i < argc; i++) {
	arg = argv[i];
	if ('-' == arg[0]) { /* it is an option */
	    if (NULL != strchr(arg,'f')) force = true;
	    if (NULL != strchr(arg,'v')) verb = true;
	    if (NULL != strchr(arg,'i')) inquire = true;
	} else { /* it is a file */
	    sf_rm(arg, force, verb, inquire);
	}
    }

    exit (0);
}
