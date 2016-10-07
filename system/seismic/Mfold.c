/* Make a seismic foldplot/stacking chart. 

This is a general 3D histogram program implemented to create foldplot or
stacking charts on a 3d project from trace headers.  Axis1, 2 and 3 
define the bins for the output fold map.  These are usually 
(offset,xline,iline), but you might want to compute some other
histogram.  This can be done by selecting other segy headers using 
label1, 2 and 3.

See also fold= option in sfbin for creating 2D histograms.

EXAMPLES:

   To make a stacking chart movie showing fold(xline,offset) for each 
   iline from a 3D segyfile:

   sfsegyread tfile=tteapot.rsf hfile=teapot.asc bfile=teapot.bin \\
           tape=npr3_field.sgy > teapot.rsf

   # read the tfile, which contains the segy trace headers
   < teapot_hdr.rsf sffold verbose=1        \\
            o1=0 n1=96  d1=200 label1=offset \\
            o2=1 n2=188 d2=1   label2=xline  \\
            o3=1 n3=345 d3=1   label3=iline  \\
   >foldplot.rsf
   <foldplot.rsf sfgrey title=foldplot pclip=100 \\
   | sfpen 

  # transpose this data to plot foldmaps for each offset window:

  < foldplot.rsf sftransp plane=13          \\
  | sftransp plane=12                       \\
  | sfgrey title=foldplot_off gainpanel=all \\
  | sfpen
  
*/
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include "segy.h"

int main(int argc, char* argv[])
{
    int verbose;
    float o1, o2, o3, d1, d2, d3;
    int n1, n2, n3, n1_input, n_input, i_input;
    int dim_in;
    int iaxis;
    char parameter[13];
    off_t n_in[SF_MAX_DIM];
    char* data_format=NULL;
    sf_datatype typehead;
    int idx_offset,idx_xline,idx_iline;
    char* label1;
    char* label2;
    char* label3;
    float* hdrin;
    float*** foldplot;
    
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");
    sf_putint(out,"input",2);

    if(!sf_getint("verbose",&verbose))verbose=1;
    /* 0 terse, 1 informative, 2 chatty, 3 debug */

    /* get o1, n1, d1, label1 for axis1. Do the same for axis2 and 3. */

    if(!sf_getfloat("o1",&o1))sf_error("o1 is a required parameter");
    /* Minimum label1 - usually min offset */

    if(!sf_getfloat("o2",&o2))sf_error("o2 is a required parameter");
    /* Minimum label2 - usually min xline  */

    if(!sf_getfloat("o3",&o3))sf_error("o3 is a required parameter");
    /* Minimum label3 - usually min iline */

    if(!sf_getint("n1",&n1))sf_error("n1 is a required parameter");
    /* Number label1 - usually number offset */

    if(!sf_getint("n2",&n2))sf_error("n2 is a required parameter");
    /* Number label2 - usually number xline */

    if(!sf_getint("n3",&n3))sf_error("n3 is a required parameter");
    /* Number label3 - usually number iline */

    if(!sf_getfloat("d1",&d1))sf_error("d1 is a required parameter");
    /* Delta label1 - usually delta offset  */

    if(!sf_getfloat("d2",&d2))sf_error("d2 is a required parameter");
    /* Delta label2 - usually delta xline  */

    if(!sf_getfloat("d3",&d3))sf_error("d3 is a required parameter");
    /* Delta label3 - usually delta iline  */

    if(NULL == (label1 = sf_getstring("label1")))label1="offset";
    /* header for axis1 - usually offset */

    if(NULL == (label2 = sf_getstring("label2")))label2="cdp";
    /*  header for axis2 - usually xline or cdp  */

    if(NULL == (label3 = sf_getstring("label3")))label3="iline";
    /* header for axis3 - usually iline  */

    /* Find out the length of headers (vectors) in the input file and the 
       number of locations to loop over 
    */

    if (!sf_histint(in,"n1",&n1_input)) 
	sf_error("input file does not define n1");

    n_input = sf_leftsize(in,1); /* left dimensions after the first one */

    data_format=sf_histstring(in,"data_format");
    fprintf(stderr,"data_format=%s\n",data_format);
    if(strcmp (data_format,"native_int")==0) typehead=SF_INT;
    else                                       typehead=SF_FLOAT;

    /* allocate space for the one location from the input file */ 
    hdrin = sf_floatalloc(n1_input);
    segy_init(n1_input,in);

    /* The output file will have new shape.  The axis have new lengths and
       names.  Write this information to the output history file.
    */
    if(verbose >=1)fprintf(stderr,"write axis information to output history\n");
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);
    sf_putfloat(out,"o3",o3);
    
    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",n3);

    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"d3",d3);
    /* if input file is more than 3 dimension, force dimensions lengths
       after 3 to be one (i.e. add n4=1 n5=1 upto iput dimension) */
    dim_in=sf_largefiledims(in,n_in);
    for (iaxis=3; iaxis<dim_in; iaxis++){
      sprintf(parameter,"n%d",iaxis+1);
      sf_putint(out,parameter,1);
    }
 
    sf_putstring(out,"label1",label1);
    sf_putstring(out,"label2",label2);
    sf_putstring(out,"label3",label3);

    /* make output file float, even it input headers are integer */
    sf_setformat(out,"native_float");

    /* get the location of the three labels in the header */
    idx_offset=segykey (label1);
    idx_xline =segykey (label2);
    idx_iline =segykey (label3);

    /* allocate the output data.  This is one big array in memory */
    foldplot=sf_floatalloc3(n1,n2,n3);

    /* loop over the input locations.  Read a trace header and increment
       the fold in each bin.  This is a 3d histogram.
    */

    if(verbose >=1)fprintf(stderr,"loop processing input trace headers\n");
    for (i_input=0 ; i_input<n_input; i_input++) {
	int iiline,ixline,ioffset;
	sf_floatread(hdrin,n1_input,in);
	/* if input are integers, make hdrin integer to read correctly */ 
	if(SF_INT == typehead){ 
	  if(verbose>2)fprintf(stderr,"typehead=SF_INT\n");
	  ioffset=roundf((((int*)hdrin)[idx_offset]-o1)/d1);
	  ixline =roundf((((int*)hdrin)[idx_xline ]-o2)/d2);
	  iiline =roundf((((int*)hdrin)[idx_iline ]-o3)/d3);
	}else{
	  if(verbose>2)fprintf(stderr,"typehead=SF_FLOAT\n");
	  ioffset=roundf((         hdrin [idx_offset]-o1)/d1);
	  ixline =roundf((         hdrin [idx_xline ]-o2)/d2);
	  iiline =roundf((         hdrin [idx_iline ]-o3)/d3);
	}

	if(verbose>2){
	    fprintf(stderr,"offset=%f,xline=%f,iline=%f\n",
		    hdrin[idx_offset],
		    hdrin[idx_xline ],
		    hdrin[idx_iline ]);
	    fprintf(stderr,"ioffset=%d,ixline=%d,iiline=%d\n",
		    ioffset   ,ixline   ,iiline);
	}
	if(ioffset>=0 && ioffset< n1 &&
	   ixline >=0 && ixline < n2 &&
	   iiline >=0 && iiline < n3 ){
	    if(verbose>2)
		fprintf(stderr,"increment fold at %d,%d,%d\n",
			iiline,ixline,ioffset); 
	    foldplot[iiline][ixline][ioffset]++;
	}
    }
    if(verbose >=1)fprintf(stderr,"write foldplot to output\n");
    /* write the output in one big write */
    sf_floatwrite(&(foldplot[0][0][0]),n1*n2*n3,out);
    
    exit(0);
}
