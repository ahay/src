/* Construct incremental minimum or maximum lists from an RSF file.

Constructs the following set of minimum or maximum lists for each
x2, x3, ... xn in the input RSF file:

out[0] = in[0]
out[i] = min or max of (in[i], out[i-1]) for i = 1, 2, 3, ... n1

The input file data type must be float.
The output file data type will be float.

sflistminmax mode=min, can be used to simulate "erosion" for a set of 
geological surfaces, producing a new set of surfaces that do not cross.

See also: sfminmax, sfstack.
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
    char        *mode=NULL;             /* Standard types           */
    int         i1,i2,n1,n2;
    float       *list=NULL;

    sf_file     in=NULL,out=NULL;       /* RSF types                */

    sf_init(argc,argv);                 /* Initialize RSF           */

    in  = sf_input("in");               /* Open files               */
    out = sf_output("out");

    sf_fileflush(out,in);               /* Copy header to output    */

                                        /* Check input data type    */
    if (sf_gettype(in) != SF_FLOAT) sf_error("Input must be float.");

                                        /* Get the file sizes       */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");

    n2 = sf_leftsize(in,1);

                                        /* Check mode               */
    if ((mode = sf_getstring("mode")) == NULL) mode = "min";
    /* 'min' (default) or 'max' */

    if (strcmp(mode,"min")!=0 && strcmp(mode,"max")!=0)
        sf_error ("Unknown mode %s.",mode);

    list = sf_floatalloc(n1);           /* Allocate list            */

    for (i2=0; i2<n2; i2++)             /* Process each list        */
    {
        sf_floatread(list,n1,in);       /* Read a list              */

                                        /* Find the min or max      */
        if (strcmp(mode,"min")==0)
            for (i1=1; i1<n1; i1++) list[i1] = SF_MIN(list[i1],list[i1-1]);

        if (strcmp(mode,"max")==0)
            for (i1=1; i1<n1; i1++) list[i1] = SF_MAX(list[i1],list[i1-1]);

        sf_floatwrite(list,n1,out);     /* Write the output         */
    }



    exit (0);
}

/* 	$Id: Mlistminmax.c 4796 2009-09-29 19:39:07Z ivlad $	 */
