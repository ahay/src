/* Element by element minimum or maximum of two RSF files.

file1 and file2 must have the same number of elements.

See also: sflistminmax, sfstack.
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

#include <rsf.h>                        /* RSF header               */
#include <string.h>

int main(int argc, char* argv[])
{
    char        *mode;                  /* Standard types           */
    long        i,n,nbuf,nread;         
    float       *v1,*v2,*v3;
    
    sf_file     file1,file2,out;        /* RSF types                */
    
    sf_init(argc,argv);                 /* Initialize RSF           */
    
                                        /* Check file parameters    */
    if (sf_getstring("file1") == NULL) sf_error("Need file1=");
    /* RSF filename required, data type must be float */
    if (sf_getstring("file2") == NULL) sf_error("Need file2=");
    /* RSF filename required, data type must be float */

    file1 = sf_input("file1");          /* Open files               */
    file2 = sf_input("file2");
    out   = sf_output("out");
    
    sf_fileflush(out,file1);            /* Copy header to output    */

                                        /* Check input data types   */
    if (sf_gettype(file1) != SF_FLOAT) sf_error("file1 must be float.");
    if (sf_gettype(file2) != SF_FLOAT) sf_error("file2 must be float.");
    
    n = sf_filesize(file1);             /* Check the file sizes     */
    
    if (sf_filesize(file2) != n)
        sf_error ("file1 and file2 must be the same size.");
        
                                        /* Check mode               */
    if ((mode = sf_getstring("mode")) == NULL) mode = "min";
    /* 'min' (default) or 'max' */ 

    if (strcmp(mode,"min")!=0 && strcmp(mode,"max")!=0)
        sf_error ("Unknown mode %s.",mode);
    
    nbuf = BUFSIZ/sizeof(float);        /* Get buffer size          */
    
    v1 = sf_floatalloc(nbuf);           /* Allocate arrays          */
    v2 = sf_floatalloc(nbuf);
    v3 = sf_floatalloc(nbuf);
    
                                        /* Read files with buffers  */
    for (nread = n; nread > 0; nread -= nbuf)
    {
        if (nbuf > nread) nbuf = nread;
        
        sf_floatread(v1,nbuf,file1);    /* Read the files           */
        sf_floatread(v2,nbuf,file2);
    
                                        /* Find the min or max      */
        if (strcmp(mode,"min")==0)
            for (i=0; i<nbuf; i++) v3[i] = SF_MIN(v1[i],v2[i]);
    
        if (strcmp(mode,"max")==0)
            for (i=0; i<nbuf; i++) v3[i] = SF_MAX(v1[i],v2[i]);
    
        sf_floatwrite(v3,nbuf,out);     /* Write the output         */
    }
    

    exit(0);
}

/* 	$Id$	 */
