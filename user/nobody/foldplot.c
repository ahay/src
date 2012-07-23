/* 
   build a fold map
   Program change history:
   date       Who             What
   12/10/2011 Karl Schleicher Original program
*/

#include <string.h>
#include <rsf.h>

#include "/home/yihua/RSFSRC/build/system/seismic/segy.h"
#include "/home/yihua/RSFSRC/build/system/seismic/segy.c"

int main(int argc, char* argv[])
{
  int verbose;
    float o1,o2,o3;
    int n1,n2,n3;
    float d1,d2,d3;
    char* label1;
    char* label2;
    char* label3;

    int n1_input;
    int n_input;
    float* hdrin;
    int idx_offset,idx_xline,idx_iline;
    int i_input;
    float*** foldplot;
    
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    sf_putint(out,"input",2);

    /* verbose flag controls ammount of print */
    /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
    if(!sf_getint("verbose",&verbose))verbose=1;

    /* Axis1, 2 and 3 define the bins for the output fold map.  These 
       are usually (offset,xline,offset), but you might want to compute some 
       other histogram.  This can be done by selecting other segy headers 
       using label1, 2 and 3.
    */

    /* get o1, n1, d1, label1 for axis1.  Do the same for axis2 and 3. */

    /*( o1=-32.5 required parameter.  Minimum label1 - usually min offset ) */
    if(!sf_getfloat("o1",&o1))sf_error("o1 is a required parameter");
    /*( o2=32 required parameter.  Minimum label2 - usually min xline ) */
    if(!sf_getfloat("o2",&o2))sf_error("o2 is a required parameter");
    /*( o3=32 required parameter.  Minimum label3 - usually min iline ) */
    if(!sf_getfloat("o3",&o3))sf_error("o3 is a required parameter");

    /*( n1=96 required parameter.  Number label1 - usually number offset ) */
    if(!sf_getint("n1",&n1))sf_error("n1 is a required parameter");
    /*( n1=623 required parameter.  Number label2 - usually number xline ) */
    if(!sf_getint("n2",&n2))sf_error("n2 is a required parameter");
    /*( n1=623 required parameter.  Number label3 - usually number iline ) */
    if(!sf_getint("n3",&n3))sf_error("n3 is a required parameter");

    /*( d1=110 required parameter.  delta label1 - usually delta offset ) */
    if(!sf_getfloat("d1",&d1))sf_error("d1 is a required parameter");
    /*( d2=1 required parameter.  delta label2 - usually delta xline ) */
    if(!sf_getfloat("d2",&d2))sf_error("d2 is a required parameter");
    /*( d3=1 required parameter.  delta label3 - usually delta iline ) */
    if(!sf_getfloat("d3",&d3))sf_error("d3 is a required parameter");

    /*( label1=offset  header for axis1 - usually offset ) */
    if(NULL == (label1 = sf_getstring("label1")))label1="offset";
    /*( label2=xline  header for axis2 - usually xline or cdp ) */
    if(NULL == (label2 = sf_getstring("label2")))label2="cdp";
    /*( label3=iline  header for axis3 - usually iline ) */
    if(NULL == (label3 = sf_getstring("label3")))label3="iline";

    if (SF_FLOAT != sf_gettype (in)) sf_error("Need float input");

    /* Find out the length of headers (vectors) in the input file and the 
       number of locations to loop over 
    */

    if (!sf_histint(in,"n1",&n1_input)) 
      sf_error("input file does not define n1");

    n_input = sf_leftsize(in,1); /* left dimensions after the first two */

    /* allocate space for the one location from the input file */ 
    hdrin = sf_floatalloc(n1_input);

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

    sf_putstring(out,"label1",label1);
    sf_putstring(out,"label2",label2);
    sf_putstring(out,"label3",label3);

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
      ioffset=roundf((hdrin[idx_offset]-o1)/d1);
      ixline =roundf((hdrin[idx_xline ]-o2)/d2);
      iiline =roundf((hdrin[idx_iline ]-o3)/d3);
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
