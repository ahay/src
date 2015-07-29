/* Trace And Header MAKEEVENT makes constant velocity dipping event synthetic.

tah is the abbreviation of Trace And Header.  Madagascar programs 
that begin with sftah are a designed to:
1- read trace and headers from separate rsf files and write them to 
   standard output (ie sftahread)
2- filter programs that read and write standard input/output and 
   process the tah data (eg sftahnmo, sftahstack)
3- read tah data from standard input and write separate rsf files for 
   the trace and headers data (ie sftahwrite)

These programs allow Seismic Unix (su) like processing in Madagascar.  
Some programs have su like names.

Some programs in this suite are sftahread, sftahgethw, ftahhdrmath, 
and sftahwrite.

The sftahmakeevent program makes simple synthetic on input data source
and group xy coordinates (i.e. sx, sy, gx, gy).  The event has constant
velocity and dip.  The nmo velocity will depend on source/receiver azimuth.
EXAMPLE:

sftahsort          \\
   input=npr3_gathers \\
   sort="iline:169,181  xline:104,116 offset:0,11000"  \\ 
   verbose=1       \\
| sftahmakeevent   \\
| sftahmakeskey pkey=iline,xline skey=cdpt verbose=1  \\      
| sftahnmo         \\
  verbose=1        \\
  tnmo=0.000,4.000 \\
  vnmo=11000,11000 \\
| sftahwrite       \\
  verbose=1        \\      
  label2="cdpt"  o2=1 n2=34  d2=1     \\
  label3="xline" o3=104 n3=13 d3=1    \\
  label4="iline" o4=169 n4=13  d4=1   \\
  output=gather_subset.rsf            \\
>/dev/null

The headers are in the file npr3_gathers_hdr.rsf, 
the headers parameter default.  The headers are merged with the trace 
amplitudes and the tah data sent down the pipe for sftahmakeevent.  The 
Traces are sent to STDOUT to have nmo applied and data to be written to 
disk.  Data is also sent to /dev/null (the bit bucket).

PARAMETERS
   none yest, but eventually ude this template:
   strings key= no default

        list of header keys print.  Look at the sfsegyread for a list
	of header names.

*/

/*
   Program change history:
   date       Who             What
   06/25/2014 Karl Schleicher Original program
*/
#include <string.h>
#include <rsf.h>
#include <rsfsegy.h>

#include "tahsub.h"

int main(int argc, char* argv[])
{
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_traces;
  int n1_headers;

  char* header_format=NULL;
  sf_datatype typehead;
  /* kls do I need to add this?  sf_datatype typein; */
  float* fheader=NULL;
  float* intrace=NULL;
  float    sx;
  float    sy;
  float    gx;
  float    gy;
  float    offset;
  float offx, offy, midx, midy;
  float scalco, scale;
  int indx_sx;
  int indx_sy;
  int indx_gx;
  int indx_gy;
  int indx_offset;
  int indx_scalco;

  float v, dx, dy, x0, y0, t0, t0xy, t, off_dotprod_dip;
  int it;
  float d1, o1;
  
 /*****************************/
  /* initialize verbose switch */
  /*****************************/
  sf_init (argc,argv);

  if(!sf_getint("verbose",&verbose))verbose=1;
      /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  sf_warning("verbose=%d",verbose);
 
  if(!sf_getfloat("v",&v))sf_error("v, a required parameter not input");
  if(!sf_getfloat("dx",&dx))sf_error("dx, a required parameter not input");
  if(!sf_getfloat("dy",&dy))sf_error("dy, a required parameter not input");
  if(!sf_getfloat("x0",&x0))sf_error("x0, a required parameter not input");
  if(!sf_getfloat("y0",&y0))sf_error("y0, a required parameter not input");
  if(!sf_getfloat("t0",&t0))sf_error("t0, a required parameter not input");

  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/

  if(verbose>0)fprintf(stderr,"read in file name\n");  
  in = sf_input ("in");

  if(verbose>0)fprintf(stderr,"read out file name\n");
  out = sf_output ("out");

  if (!sf_histint(in,"n1_traces",&n1_traces))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n1_headers",&n1_headers)) 
    sf_error("input data does not define n1_headers");

  header_format=sf_histstring(in,"header_format");
  if(strcmp (header_format,"native_int")==0) typehead=SF_INT;
  else                                       typehead=SF_FLOAT;

  if(verbose>0)fprintf(stderr,"allocate headers.  n1_headers=%d\n",n1_headers);
  fheader = sf_floatalloc(n1_headers);
 
  if(verbose>0)fprintf(stderr,"allocate intrace.  n1_traces=%d\n",n1_traces);
  intrace= sf_floatalloc(n1_traces);



  /* maybe I should add some validation that n1== n1_traces+n1_headers+2
     and the record length read in the second word is consistent with 
     n1.  */

  /**********************************************************/
  /* end code block for standard tah Trace And Header setup */
  /* continue with any sf_puthist this tah program calls to */
  /* add to the history file                                */
  /**********************************************************/

  /* put the history from the input file to the output */
  sf_fileflush(out,in);

  /********************************************************/
  /* continue initialization specific to this tah program */
  /********************************************************/
  if (!sf_histfloat(in,"d1",&d1))
    sf_error("input data not define d1");
  if (!sf_histfloat(in,"o1",&o1))
    sf_error("input data not define o1");

  /* segy_init gets the list header keys required by segykey function  */
  segy_init(n1_headers,in);

  indx_sx    =segykey("sx"    );
  indx_sy    =segykey("sy"    );
  indx_gx    =segykey("gx"    );
  indx_gy    =segykey("gy"    );
  indx_offset=segykey("offset");
  indx_scalco=segykey("scalco");

  /***************************/
  /* start trace loop        */
  /***************************/
  if(verbose>0)fprintf(stderr,"start trace loop\n");
  while (0==get_tah(intrace, fheader, n1_traces, n1_headers, in)){
    if(verbose>1)fprintf(stderr,"process the tah in sftahgethw\n");
    /********************/
    /* process the tah. */
    /********************/
    if(typehead == SF_INT){
      /* just cast the header to int so the print works */
      scalco=((int*)fheader)[indx_scalco];
      if(scalco==0)scale=1;
      if(scalco<0)scale=-1/scalco;
      else scale=scalco;
      sx    =scale*((int*)fheader)[indx_sx    ];
      sy    =scale*((int*)fheader)[indx_sy    ];
      gx    =scale*((int*)fheader)[indx_gx    ];
      gy    =scale*((int*)fheader)[indx_gy    ];
      offset=((int*)fheader)[indx_offset];
    } else {
      scalco=fheader[indx_scalco];
      if(scalco==0)scale=1;
      if(scalco<0)scale=-1/scalco;
      else scale=scalco;
      sx    =scale*fheader[indx_sx    ];
      sy    =scale*fheader[indx_sy    ];
      gx    =scale*fheader[indx_gx    ];
      gy    =scale*fheader[indx_gy    ];
      offset=fheader[indx_offset];
    }
    offx=sx-gx;
    offy=sy-gy;
    midx=(sx+gx)/2.0;
    midy=(sy+gy)/2.0;
    if(verbose>2){
      fprintf(stderr, "sx    =%e",sx    );
      fprintf(stderr," sy    =%e",sy    );
      fprintf(stderr," gx    =%e",gx    );
      fprintf(stderr," gy    =%e",gy    );
      fprintf(stderr," offset=%e",offset);
      fprintf(stderr," offx=%e",offx);
      fprintf(stderr," offy=%e",offy);
      fprintf(stderr," newoff=%e",sqrt(offx*offx+offy*offy));
      fprintf(stderr," midx=%e",midx);
      fprintf(stderr," midy=%e",midy);
      fprintf(stderr,"\n");
    }
    /* zero the output trace */
    memset(intrace,0,n1_traces*sizeof(float));
    /* compute time of the spike and put it on the trace */
    /* t=t0; */
    t0xy=t0+(midx-x0)*dx+(midy-y0)*dy;
    /* t=t0xy; */
    /* no dmo term t=sqrt(t0xy*t0xy+(offset*offset)/(v*v)); */
    off_dotprod_dip=offx*dx+offy*dy;
    t=sqrt(t0xy*t0xy+(offset*offset)/(v*v)-off_dotprod_dip*off_dotprod_dip/4.0);
    if(t>=o1){
      it=(t-o1)/d1;
      if(it>0 && it<n1_traces){
	intrace[it]=it+1-(t-o1)/d1;
      }
      if(it+1>0 && it+1<n1_traces){
	intrace[it+1]=(t-o1)/d1-it;
      }
    }
    /***************************/
    /* write trace and headers */
    /***************************/
    put_tah(intrace, fheader, n1_traces, n1_headers, out);
  }

  exit(0);
}

  
