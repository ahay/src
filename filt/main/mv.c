#include <string.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i;
    char *infile = NULL;
    sf_file in = NULL, out = NULL;

    sf_init (argc,argv);

    /* the first two non-parameters are in and out files */
    for (i=1; i< argc; i++) { 
	if (NULL == strchr(argv[i],'=')) {
	    if (NULL == in) {
		infile = argv[i];
		in = sf_input (infile);
	    } else {
		out = sf_output (argv[i]);
		break;
	    }
	}
    }
    if (NULL == in || NULL == out)
	sf_error ("not enough input");

    sf_setformat(out,sf_histstring(in,"data_format"));

    sf_cp(in,out);
    if (NULL != strstr (sf_getprog(),"mv")) sf_rm(infile,false,false,false);
    exit (0);
}
