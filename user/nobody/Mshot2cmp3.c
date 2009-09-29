/* Convert shots to CMPs for regular 3-D geometry. 

The axes in the input are {time, offset_x, offset_y, shot_x, shoty}
The axes in the output are {time, offset_x, offset_y, midpoint_x, midpoint_y}
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

int main(int argc, char* argv[])
{

    off_t pos;
    bool sign, half;
    int   nsx, nsy,   nmx, nmy,   nhx, nhy,   nmhx, nmhy,   nt;
    int   isx, isy,   imx, imy,   ihx, ihy;
    float osx, osy, dsx, dsy, omx, omy, dmx, dmy, ohx, ohy, dhx, dhy, dmhx, dmhy, binx, biny;
    float s_xmin, s_ymin, r_xmin, r_ymin, s_xmax, s_ymax, r_xmax, r_ymax, survey_xmin, survey_ymin, survey_xmax, survey_ymax;
    char *trace, *zero;
    sf_file in, out;


    sf_init(argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");
    sf_warning("WARNING: This 3-D CMP sorting code is not yet ready for use--this message will be removed when it is. ");

    if (!sf_histint  (in,"n1",&nt)) sf_error("No n1= in input");
    /* Number of samples per trace */

    if (!sf_histint  (in,"n2",&nhx)) sf_error("No n2= in input");
    /* Number of offsets per shot gather in x-direction*/
    if (!sf_histfloat(in,"o2",&ohx)) sf_error("No o2= in input");
    /* First offset x-component */
    if (!sf_histfloat(in,"d2",&dhx)) sf_error("No d2= in input");
    /* Offset increment in x-direction */
    if (!sf_histint  (in,"n3",&nhy)) sf_error("No n3= in input");
    /* Number of offsets per shot gather in y-direction*/
    if (!sf_histfloat(in,"o3",&ohy)) sf_error("No o3= in input");
    /* First offset y-component */
    if (!sf_histfloat(in,"d3",&dhy)) sf_error("No d3= in input");
    /* Offset increment in y-direction */
    if (!sf_histint  (in,"n4",&nsx)) sf_error("No n4= in input");
    /* Number of sources along x-direction*/
    if (!sf_histfloat(in,"d4",&dsx)) sf_error("No d4= in input");
    /* Source spacing in x-direction*/
    if (!sf_histfloat(in,"o4",&osx)) sf_error("No o4= in input");
    /* First source x-coordinate*/
    if (!sf_histint  (in,"n5",&nsy)) sf_error("No n5= in input");
    /* Number of sources along y-direction*/
    if (!sf_histfloat(in,"d5",&dsy)) sf_error("No d5= in input");
    /* Source spacing in y-direction*/
    if (!sf_histfloat(in,"o5",&osy)) sf_error("No o5= in input");
    /* First source y-coordinate*/
    
    if (!sf_getfloat("binx",&binx)) sf_error("No x bin size specified. Please set binx=");
    /*Number of bins along x-direction*/
    if (!sf_getfloat("biny",&biny)) sf_error("No y bin size specified. Please set biny=");
    /*Number of bins along y-direction*/

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation:
       yes is generally for off-end surveys, where the first offsets are positive.  
       no is generally for split-spread surveys with first negative then positive offsets. */

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset*/

    if (!half) {
	dhx /= 2;
	ohx /= 2;
	dhy /= 2;
	ohy /= 2;
    }

    s_xmin = osx;
    s_ymin = osy;
    s_xmax = osx+(nsx*dsx);
    s_ymax = osy+(nsy*dsy);
    
    r_xmin = s_xmin+2*ohx;
    r_ymin = s_ymin+2*ohy;
    r_xmax = s_xmax+2*(ohx+dhx*nhx);
    r_ymax = s_ymax+2*(ohy+dhy*nhy);

    if (s_xmin <= r_xmin){
      survey_xmin = s_xmin;
    }else{
      survey_xmin = r_xmin;
    }

    if (s_ymin <= r_ymin){
      survey_ymin = s_ymin;
    }else{
      survey_ymin = r_ymin;
    }

    if (s_xmax <= r_xmax){
      survey_xmax = r_xmax;
    }else{
      survey_xmax = s_xmax;
    }

    if (s_ymax <= r_ymax){
      survey_ymax = r_ymax;
    }else{
      survey_ymax = s_ymax;
    }
    
    dmx = (survey_xmax - survey_xmin)/binx;
    dmy = (survey_ymax - survey_ymin)/biny;
    omx = survey_xmin + dmx/2.;
    omy = survey_ymin + dmy/2.;
    nmx = binx;
    nmy = biny;
    nmhx = (int)(nhx*nsx/nmx);
    nmhy = (int)(nhy*nsy/nmy);
    dmhx = dhx;
    dmhy = dhy;

    sf_putint  (out,"n2",nmhx);
    sf_putint  (out,"n3",nmhy);
    sf_putfloat(out,"d2",dmhx);
    sf_putfloat(out,"d3",dmhy);

    sf_putint  (out,"n4",nmx);
    sf_putint  (out,"n5",nmy);
    sf_putfloat(out,"d4",dmx);
    sf_putfloat(out,"d5",dmy);
    sf_putfloat(out,"o4",omx);
    sf_putfloat(out,"o5",omy);

    sf_putstring(out,"label3","hx");
    sf_putstring(out,"label4","hy");
    sf_putstring(out,"label5","mx");
    sf_putstring(out,"label6","my");
    
    nt *= sf_esize(in);

    trace = sf_charalloc(nt);
    zero  = sf_charalloc(nt);
    memset(zero,0,nt);

    sf_fileflush(out,in);
    sf_setform( in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);
    
    sf_unpipe(in,(off_t) nsx*nsy*nhx*nhy*nt);
    pos = sf_tell(in);

    for (isx=0; isx < nsx; isx++) {
      for (isy=0; isy < nsy; isy++) {

	for (ihx=0; ihx < nhx; ihx++) {
	  for (ihy=0; ihy < nhy; ihy++) {

	    imx = (int)((isx + ihx)/2.0);
	    imy = (int)((isy + ihy)/2.0);

	    if (isx >= 0 && isx < nsx && ihx < nhx) {
	      if (isy >= 0 && isy < nsy && ihy < nhy) {
		sf_seek(in,pos+((isx*nhx+ihx)+(isy*nhy+ihy))*nt,SEEK_SET);
		sf_charread(trace,nt,in);
		sf_charwrite(trace,nt,out);
	      }
	    } else {
		sf_charwrite(zero,nt,out);
	    }
	  }
	}
      }
    }


    exit(0);
}
