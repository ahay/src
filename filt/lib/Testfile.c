#include <stdlib.h>

#include "file.h"
#include "getpar.h"

int main(int argc, char* argv[])
{
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    sf_fileclose (in);
    sf_fileclose (out);

    exit (0);
}

/* 	$Id: Testfile.c,v 1.2 2003/09/29 14:34:55 fomels Exp $	 */

