/* Display basic information about RSF files.

Takes: [<file0.rsf] file1.rsf file2.rsf ...

n1,n2,... are data dimensions
o1,o2,... are axis origins
d1,d2,... are axis sampling intervals
label1,label2,... are axis labels
unit1,unit2,... are axis units
*/
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
#include <string.h>

#include <rsf.h>

static void check_zeros (sf_file file, int esize, off_t size, 
			 int ncheck, char buf[], char zero[]);

int main (int argc, char* argv[])
{
    int i, j, j1, ncheck, esize, nin, dim=SF_MAX_DIM;
    off_t nj, size, n[SF_MAX_DIM];
    float check, fj;
    char *label, *dataname, key[7], *val;
    const char **filename;
    char buf[BUFSIZ], zero[BUFSIZ];
    sf_file file;
    bool info, trail;
    const char *type[] = {"uchar","char","int","float","complex","short","double","long"};
    const char *form[] = {"ascii","xdr","native"};
    char pad[] = "              ", out[30];

    sf_init (argc,argv);

    filename = (const char**) sf_alloc ((size_t) argc,sizeof(char*));

    if (!sf_getbool ("info",&info)) info = true;
    /* If n, only display the name of the data file. */
    if (!sf_getfloat ("check",&check)) check = 2.;
    /* Portion of the data (in Mb) to check for zero values. */
    check *= (1024. * 1024.); /* convert Mb to b */
    ncheck = (int) check;

    if (!sf_getbool("trail",&trail)) trail=true;
    /* If n, skip trailing dimensions of  one */

    memset(zero,0,BUFSIZ);

    if (!sf_stdin()) { /* no input file in stdin */
	nin=0;
    } else {
	filename[0] = "in";
	nin=1;
    }

    for (i = 1; i < argc; i++) {
	if (NULL != strchr (argv[i], '=')) 
	  continue; /* not a file */

	filename[nin] = argv[i];
	nin++;
    }
    if (0==nin) sf_error ("no input");
	
    for (i = 0; i < nin; i++) {
	file = sf_input (filename[i]);
	dataname = sf_histstring(file,"in");

	if (!info) {
	    printf("%s ",dataname);
            sf_fileclose(file);
	    continue;
	}

	printf ("%s:\n", filename[i]);
	printf ("%sin=\"%s\"\n",pad+10,dataname);

	if (sf_histint(file,"esize",&esize)) {
	    printf ("%sesize=%d ",pad+10,esize);
	} else {
	    esize = sf_esize(file);
	    printf ("%sesize=%d? ",pad+10,esize);
	}
	printf("type=%s form=%s ",
	       type[sf_gettype(file)],
	       form[sf_getform(file)]);

	if (NULL != (label = sf_histstring(file,"label"))) {
	    printf("label=\"%s\" ",label);
	    free(label);
	}

	if (NULL != (label = sf_histstring(file,"unit"))) {
	    printf("unit=\"%s\"",label);
	    free(label);
	}

	printf("\n");

	if (!trail) dim = sf_largefiledims(file,n);

	size = 1;
	for (j=0; j < dim; j++) {
	    j1 = (j+1)%10;
	    
	    snprintf(key,3,"n%d",j1);
	    if (!sf_histlargeint(file,key,&nj)) break;

#if defined(__cplusplus) || defined(c_plusplus)
	    snprintf(out,28,"%s=%lu",key,(long) nj);
#else
	    snprintf(out,28,"%s=%llu",key,(long long) nj);
#endif
	    printf("%s%s%s",pad+10,out,pad+strlen(out));
	    size *= nj;

	    snprintf(key,3,"d%d",j1);
	    if (sf_histfloat(file,key,&fj)) {
		snprintf(out,28,"%s=%g ",key,fj);
	    } else {
		snprintf(out,28,"%s=? ",key);
	    }
	    printf(" %s%s",out,pad+strlen(out));

	    snprintf(key,3,"o%d",j1);
	    if (sf_histfloat(file,key,&fj)) {
		snprintf(out,28,"%s=%g ",key,fj);
	    } else {
		snprintf(out,28,"%s=? ",key);
	    }
	    printf(" %s%s",out,pad+strlen(out));

	    snprintf(key,7,"label%d",j1);
	    if (NULL != (val = sf_histstring(file,key))) {
		printf("%s=\"%s\" ",key,val);
	    }

	    snprintf(key,6,"unit%d",j1);
	    if (NULL != (val = sf_histstring(file,key))) {
		printf("%s=\"%s\" ",key,val);
	    }

	    printf("\n");
	}

	check_zeros (file, esize, size, ncheck, buf, zero);
	sf_fileclose(file);
    }
    if (!info) printf("\n");

    exit (0);
}

static void check_zeros (sf_file file, int esize, off_t size, int ncheck,
			 char buf[], char zero[])
{
    off_t bytes;
    int nzero, nleft, nbuf;

#if defined(__cplusplus) || defined(c_plusplus)
    if (0==esize) {
	printf("\t%ld elements\n",(long int) size);
    } else {
	printf("\t%ld elements %ld bytes\n",(long int) size,
	       (long int) size * esize);
#else
    if (0==esize) {
	printf("\t%lld elements\n",(long long int) size);
    } else {
	printf("\t%lld elements %lld bytes\n",(long long int ) size,
	       (long long int) size * esize);	
#endif
	bytes = sf_bytes(file);
	size *= esize;
	
	if (bytes < 0) bytes = size;

	if (size != bytes) 
	    sf_warning("\t\tActually %lld bytes, %g%% of expected.",
		       bytes, (100.00*bytes)/size);	
	
	for (nzero=0, nleft = bytes, nbuf = BUFSIZ; 
	     nzero < ncheck && nleft > 0; 
	     nleft -= nbuf, nzero += nbuf) {
	    if (nbuf > nleft) nbuf = nleft;
	    sf_charread(buf,nbuf,file);

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
