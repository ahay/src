/* RSF file to Text (ASCII) file (as a matrix) 
Example: 
sftxt2rsf tfile=location.txt n1=2 n2=3 >location.rsf
sfrsf2txt <location.rsf format=.4f tfile=location2.txt
more location.txt
more location2.txt
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
#include <string.h>

int main(int argc, char *argv[])
{
    int i1,i2,n1,n2,n123;       		/* number of samples in input file    	*/
    float *u;	    		/* trace for storing data 	 	*/
    char *tfilename=NULL;	/* name string of output ascii file 	*/
    char *format;	/* fprintf format (default: f) 	*/
    char format2[10];	/*%f */
    FILE *tfile;		/* File pointer of output ascii file 	*/
    sf_file in;			/* input rsf file 			*/

/***************************************************/
/*	Initialization				   */
/***************************************************/
    sf_init(argc,argv);		
    in = sf_input("in");

/***************************************************/
/*	Allocating memory			   */
/***************************************************/
    
    sf_histint(in,"n1",&n1);
    n2=sf_leftsize(in,1);
    n123=n1*n2;
    u=sf_floatalloc(n123);

/***************************************************/
/*	Reading rsf file		  	    */
/***************************************************/
    sf_floatread(u,n123,in);

	strcpy(format2,"%");	
    if(NULL==(format=sf_getstring("format")))
	{
	format="f";
	}
	strcat(format2,format);
	strcat(format2," ");
	
/***************************************************/
/*	Writing ascii file		 	   */
/***************************************************/
    if(NULL==(tfilename=sf_getstring("tfile")))
	{
	tfile=stdout;	
	}
    else if(NULL==(tfile=fopen(tfilename,"wb")))
	{
	sf_error("Cannot open \"%s\" for writing:",tfilename);
	}
    
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
   			fprintf(tfile,format2,u[i1+i2*n1]);
   		if(i1<n1-1)
   		{fprintf(tfile,"\n");}
   	}
    
    fclose(tfile);

    exit(0);
}
