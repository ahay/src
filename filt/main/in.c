/* Display basic information about RSF files.
   
Takes: file1.rsf file2.rsf file3.rsf

n1,n2,... are data dimensions
o1,o2,... are axis origins
d1,d2,... are axis sampling intervals
label1,label2,... are axis labels
*/

#include <string.h>

#include <rsf.h>

static void check_zeros (sf_file file, int esize, int size, int ncheck,
			 char buf[], char zero[]);

int main (int argc, char* argv[])
{
    int i, j, ncheck, esize, dim, nj, size;
    float check, fj;
    char *filename, *dataname, key[8], *val, buf[BUFSIZ], zero[BUFSIZ];
    sf_file file;
    bool info;
    char *type[] = {"char","int","float","complex"};
    char *form[] = {"ascii","xdr","native"};

    sf_init (argc,argv);
    if (!sf_getbool ("info",&info)) info = true;
    /* If n, only display the name of the data file. */
    if (!sf_getfloat ("check",&check)) check = 2.;
    /* Portion of the data (in Mb) to check for zero values. */
    check *= (1024. * 1024.); /* convert Mb to b */
    ncheck = (int) check;

    memset(zero,0,BUFSIZ);
    
    for (i = 1; i < argc; i++) {
	filename = argv[i];
	if (NULL != strchr (filename, '=')) continue;
	
	file = sf_input (filename);
	dataname = sf_histstring(file,"in");

	if (!info) {
	    printf("%s ",dataname);
	    continue;
	}

	printf ("%s:\n", filename);
	printf ("\tin=\"%s\"\n",dataname);

	if (sf_histint(file,"esize",&esize)) {
	    printf ("\tesize=%d ",esize);
	} else {
	    esize = 4;
	    printf ("\tesize=%d? ",esize);
	}
	printf("type=%s form=%s\n",
	       type[sf_gettype(file)],
	       form[sf_getform(file)]);
	
	printf("\t");
	size = 1;
	for (j=0; j < SF_MAX_DIM; j++) {
	    sprintf(key,"n%d",j+1);
	    if (!sf_histint(file,key,&nj)) break;
	    printf("%s=%d ",key,nj);
	    size *= nj;
	}
	printf("\n");
	dim = j;

	check_zeros (file, esize, size, ncheck, buf, zero);
		
	printf("\t");
	for (j=0; j < dim; j++) {
	    sprintf(key,"d%d",j+1);
	    if (sf_histfloat(file,key,&fj)) {
		printf("%s=%g ",key,fj);
	    } else {
		printf("%s=? ",key);
	    }
	}
	printf("\n");

	printf("\t");
	for (j=0; j < dim; j++) {
	    sprintf(key,"o%d",j+1);
	    if (sf_histfloat(file,key,&fj)) {
		printf("%s=%g ",key,fj);
	    } else {
		printf("%s=? ",key);
	    }
	}
	printf("\n");

	printf("\t");
	for (j=0; j < dim; j++) {
	    sprintf(key,"label%d",j+1);
	    if (NULL != (val = sf_histstring(file,key))) {
		printf("%s=\"%s\" ",key,val);
	    } 
	}
	printf("\n");
    }
    if (!info) printf("\n");
    
    exit (0);
}

static void check_zeros (sf_file file, int esize, int size, int ncheck,
			 char buf[], char zero[])
{
    long bytes;
    int nzero, nleft, nbuf;

    if (0==esize) {
	printf("\t\t%d elements\n",size);
    } else {
	printf("\t\t%d elements %d bytes\n",size,size*esize);
	bytes = sf_bytes(file);
	size *= esize;
	
	if (bytes < 0) bytes = size;

	if (size != bytes) 
	    sf_warning("\t\tActually %d bytes, %g%% of expected.",
		       bytes, 100.*bytes/size);	
	
	for (nzero=0, nleft = bytes, nbuf = BUFSIZ; 
	     nzero < ncheck && nleft > 0; 
	     nleft -= nbuf, nzero += nbuf) {
	    if (nbuf > nleft) nbuf = nleft;
	    sf_read(buf,1,nbuf,file);
	    
	    if (0 != memcmp(buf,zero,nbuf)) break;
	}
	
	if (nzero > 0) {
	    if (nzero == size) {
		sf_warning("This data file is entirely zeros.");
	    } else if (nzero == ncheck) {
		sf_warning("This data file might be all zeros"
			   " (checked %d bytes)",nzero);
	    } else {
		sf_warning("The first %d bytes are all zeros", nzero);
	    }
	}
    }
}

/* 	$Id: in.c,v 1.4 2004/02/14 07:20:09 fomels Exp $	 */

