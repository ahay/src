/* Text (ASCII) file (like a matrix) to RSF file 
Convert a file containing a two dimensional array of ASCII floats to 
.rsf format.  n1*n2*sizeof(float) bytes are read from the input file
and the rsf file is written to standard output.  If you have a higher
dimension file (3d, 4d, ..) you can change n2, n3, etc using sfput.  

Example that converts and plots a ascii velocity grid:

<Vp sftxt2rsf n1=400 d1=1 n2=1600 d2=1 o2=101 \\
| sfgrey color=jet scalebar=y allpos=y | sfpen
 
*/
/*
  Copyright (C) 2021 University of Texas at Austin
   
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

#include <rsf.h>

int main(int argc, char *argv[])
{
	int i1,i2;
    int n1;       		/* number of samples in first dimension */
    int n2;       		/* number of samples in second dimension*/
    float d1;       		/* temporal interval in first dimension */
    float d2;       		/* temporal interval in second dimension*/
    float o1;       		/* starting position in first dimension */
    float o2;       		/* starting position in second dimension*/
    float *u;	    		/* trace for storing data 	 	*/
    char *tfilename=NULL;	/* name string of output ascii file 	*/
    FILE *tfile;		/* File pointer of output ascii file 	*/
    sf_file out;		/* input rsf file 			*/

/***************************************************/
/*	Initialization				   */
/***************************************************/
    sf_init(argc,argv);		
    out = sf_output("out");

/***************************************************/
/*	Getting and putting dimensions 		   */
/***************************************************/
    if(!sf_getint("n1",&n1)) sf_error(" No n1 specified !");
    if(!sf_getint("n2",&n2)) sf_error(" No n2 specified !");
    if(!sf_getfloat("d1",&d1)) d1=0.004;
    if(!sf_getfloat("d2",&d2)) d2=1;
    if(!sf_getfloat("o1",&o1)) o1=0;
    if(!sf_getfloat("o2",&o2)) o2=0;

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);

/***************************************************/
/*		Allocate memory 		   */
/***************************************************/
    u=sf_floatalloc(n1*n2);

/***************************************************/
/*		Reading ascii file		   */
/***************************************************/
    if(NULL==(tfilename=sf_getstring("tfile")))
	{
	tfile=stdin;	
	}
    else if(NULL==(tfile=fopen(tfilename,"rb")))
	{
	sf_error("Cannot open \"%s\" for reading:",tfilename);
	}
    /* if tfile is not input, data will be read from stdin */
   
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
   			fscanf(tfile,"%f",u+i1+i2*n1);
   	}
	
	fclose(tfile);
/***************************************************/
/*	Writing rsf file		   */
/***************************************************/
    sf_floatwrite(u,n1*n2,out);

    exit(0);
}
