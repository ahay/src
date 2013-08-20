/* RSF file to Binary file */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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
    int n123;       		/* number of samples in input file    	*/
    float *u;	    		/* trace for storing data 	 	*/
    char *bfilename=NULL;	/* name string of output binary file 	*/
    FILE *bfile;		/* File pointer of output binary file 	*/
    sf_file in;			/* input rsf file 			*/

/***************************************************/
/*	Initialization				   */
/***************************************************/
    sf_init(argc,argv);		
    in = sf_input("in");

/***************************************************/
/*	Allocating memory			   */
/***************************************************/
    n123 = sf_leftsize(in,0);
    u=sf_floatalloc(n123);

/***************************************************/
/*	Reading rsf file		  	    */
/***************************************************/
    sf_floatread(u,n123,in);

/***************************************************/
/*	Writing binary file		 	   */
/***************************************************/
    if(NULL==(bfilename=sf_getstring("bfile")))
	{
	bfile=stdout;	
	}
    else if(NULL==(bfile=fopen(bfilename,"wb")))
	{
	sf_error("Cannot open \"%s\" for writing:",bfilename);
	}
    fwrite(u,1,n123*sizeof(float),bfile);
    fclose(bfile);

    exit(0);
}
