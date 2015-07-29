#include <stdlib.h>
#include <assert.h>

#include "file.h"
#include "getpar.h"

int main(int argc, char* argv[])
{
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    assert (NULL != in);
    assert (NULL != out);

    sf_fileclose (out);

    exit (0);
}
