/* Bilinear interpolation of a 2-D array of bytes */
/*
  Copyright © 2007, Colorado School of Mines,
  All rights reserved.
  
  
  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the following 
  conditions are met:
  
  *  Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  *  Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions and the following 
  disclaimer in the documentation and/or other materials provided 
  with the distribution.
  *  Neither the name of the Colorado School of Mines nor the names of
  its contributors may be used to endorse or promote products 
  derived from this software without specific prior written permission.
  
  Warranty Disclaimer:
  THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  POSSIBILITY OF SUCH DAMAGE.
  
  
  Export Restriction Disclaimer:
  We believe that CWP/SU: Seismic Un*x is a low technology product that does
  not appear on the Department of Commerce CCL list of restricted exports.
  Accordingly, we believe that our product meets the qualifications of
  an ECCN (export control classification number) of EAR99 and we believe
  it fits the qualifications of NRR (no restrictions required), and
  is thus not subject to export restrictions of any variety.
  
  Approved Reference Format:
  In publications, please refer to SU as per the following example:
  Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
  Release No. __: an open source software  package for seismic 
  research and processing, 
  Center for Wave Phenomena, Colorado School of Mines.
  
  Articles about SU in peer-reviewed journals:
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/
#include <rsf.h>

/* number of tabulated interpolation coefficients */
#define ICMAX 99 /* must be odd, so that ICMAC-ic!=ic, for ic=0 to ICMAX/2! */
#define NTABLE (ICMAX+1)

/* functions defined and used internally */
static void intl2bx(int, int*, int*, int,
		    unsigned char[][256], unsigned char*, unsigned char*);
static void intl2by(int, int, int, unsigned char[][256],
		    unsigned char*, unsigned char*, unsigned char*);

void intl2b(int nxin            /* number of x samples input (fast dimension of zin) */, 
	    float dxin          /* x sampling interval input */, 
	    float fxin          /* first x sample input */,
	    int nyin            /* number of y samples input (slow dimension of zin) */, 
	    float dyin          /* y sampling interval input */, 
	    float fyin          /* first y sample input */, 
	    unsigned char *zin  /* array[nyin][nxin] of input samples (see notes) */,
	    int nxout           /* number of x samples output (fast dimension of zout) */, 
	    float dxout         /* x sampling interval output */, 
	    float fxout         /* first x sample output */,
	    int nyout           /* number of y samples output (slow dimension of zout) */, 
	    float dyout         /* y sampling interval output */, 
	    float fyout         /* first y sample output */, 
	    unsigned char *zout /* [nyout][nxout] output samples (see notes) */)
/*< bilinear interpolation of a 2-D array of bytes 
  Notes:
  The arrays zin and zout must passed as pointers to the first element of
  a two-dimensional contiguous array of unsigned char values.

  Constant extrapolation of zin is used to compute zout for
  output x and y outside the range of input x and y.
  
  For efficiency, this function builds a table of interpolation
  coefficents pre-multiplied by byte values.  To keep the table
  reasonably small, the interpolation does not distinguish
  between x and y values that differ by less than dxin/ICMAX
  and dyin/ICMAX, respectively, where ICMAX is a parameter
  defined above. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/02/89
Modified:  Dave Hale, Colorado School of Mines, 05/30/90
	Changed function to interpolate unsigned char 
	instead of signed char, since many color tables and 
	image processing functions (e.g., PostScript) require 
	bytes with a maximum range of 0 to 255.
Modified:  Dave Hale, Colorado School of Mines, 06/01/91
	Changed computation of coefficient table to avoid
	errors due to truncation in float to fix.  Old code
	sometimes caused interpolated values to be less than
	the minimum of the byte values being interpolated or 
	greater than the maximum of the values being interpolated.
*****************************************************************************/
{         
    int ixout,iyout,ic,ib,iyin,iyinl=1;
    float xout,yout,rxin,ryin,frac;
    int *kzin,*kic;
    unsigned char *temp1,*temp2,*temp;
    static unsigned char table[NTABLE][256];
    static int tabled=0;

    /* if not already built, build byte multiplication table */
    if (!tabled) {
	for (ic=0; ic<=ICMAX/2; ++ic) {
	    frac = (float)(ic)/(float)ICMAX;
	    for (ib=0; ib<256; ++ib) {
		table[ic][ib] = frac*ib;
		table[ICMAX-ic][ib] = ib-table[ic][ib];
	    }
	}
	tabled = 1;
    }                                              

    /* get workspace */
    kzin = sf_intalloc(nxout);
    kic = sf_intalloc(nxout);
    temp1 = sf_ucharalloc(nxout);
    temp2 = sf_ucharalloc(nxout);

    /* pre-compute indices for fast 1-D interpolation along x axis */
    for (ixout=0,xout=fxout; ixout<nxout; ixout++,xout+=dxout) {
	rxin = (xout-fxin)/dxin;
	if (rxin<=0) {
	    kzin[ixout] = 0;
	    kic[ixout] = 0;
	} else if (rxin>=nxin-1) {
	    kzin[ixout] = nxin-2;
	    kic[ixout] = ICMAX*256;
	} else {
	    kzin[ixout] = (int)rxin;
	    frac = rxin-(int)rxin;
	    ic = frac*ICMAX+0.5;
	    kic[ixout] = ic*256;
	}
    }

    /* loop over output y */
    for (iyout=0,yout=fyout; iyout<nyout; iyout++,yout+=dyout) {

	/* compute index of input y, clipped to range of input y */
	ryin = SF_MAX(0,SF_MIN(nyin-1,(yout-fyin)/dyin));
	iyin = SF_MAX(0,SF_MIN(nyin-2,ryin));

	/* if output y is not between current input y */
	if (iyin!=iyinl || iyout==0) {

	    /* if 2nd temporary vector is still useful */
	    if (iyin==iyinl+1 && iyout!=0) {              

		/* swap 2nd and 1st temp; compute 2nd temp */
		temp = temp1;
		temp1 = temp2;
		temp2 = temp;
		intl2bx(nxout,kzin,kic,ICMAX,
			table,zin+(iyin+1)*nxin,temp2);

		/* else if 1st temporary vector is still useful */
	    } else if (iyin==iyinl-1 && iyout!=0) {

		/* swap 1st and 2nd temp; compute 1st temp */
		temp = temp1;
		temp1 = temp2;
		temp2 = temp;
		intl2bx(nxout,kzin,kic,ICMAX,
			table,zin+iyin*nxin,temp1);

		/* else if neither 1st or 2nd temp is useful */
	    } else {

		/* compute 1st and 2nd temporary vectors */
		intl2bx(nxout,kzin,kic,ICMAX,
			table,zin+iyin*nxin,temp1);
		intl2bx(nxout,kzin,kic,ICMAX,
			table,zin+(iyin+1)*nxin,temp2);
	    }                 

	    /* remember last index of input y */
	    iyinl = iyin;
	}

	/* compute index of interpolation coefficient */
	frac = ryin-iyin;
	ic = frac*ICMAX+0.5;

	/* linearly interpolate output vector by table lookup */
	intl2by(nxout,ic,ICMAX,table,
		temp1,temp2,zout+iyout*nxout);
    }                         

    /* free workspace before returning */
    free(kzin);
    free(kic);                     
    free(temp1);
    free(temp2);
}
   
static void intl2bx(int nxout, int *kzin, int *kic, int icmax, 
		    unsigned char table[][256], unsigned char *zin, unsigned char *zout)
/****************************************************************************
 interpolate between input x values (FOR INTERNAL USE by intl2b) 
****************************************************************************/
{
    int ixout,jzin,jic;
    unsigned char *tablel,*tableh;
    tablel = &table[0][0];
    tableh = &table[icmax][0];
    for (ixout=0; ixout<nxout; ixout++) {   
	jzin = kzin[ixout];
	jic = kic[ixout];
	zout[ixout] = tableh[(int)zin[jzin]-jic] 
	    + tablel[(int)zin[jzin+1]+jic];
    }
}   

static void intl2by(int nxout, int ic, int icmax, unsigned char table[][256],
		    unsigned char *temp1, unsigned char *temp2, unsigned char *zout)
/****************************************************************************
 interpolate between input y values (FOR INTERNAL USE by intl2b) 
****************************************************************************/
{
    int ixout;
    unsigned char *tablel, *tableh;
    tablel = &table[ic][0];
    tableh = &table[icmax-ic][0];
    for (ixout=0; ixout<nxout; ixout++)
	zout[ixout] = tableh[temp1[ixout]]+tablel[temp2[ixout]];
}

void intl2b_block(int nxin            /* number of x samples input (fast dimension of zin) */, 
		  float dxin          /* x sampling interval input */, 
		  float fxin          /* first x sample input */,
		  int nyin            /* number of y samples input (slow dimension of zin) */, 
		  float dyin          /* y sampling interval input */, 
		  float fyin          /* first y sample input */, 
		  unsigned char *zin  /* array[nyin][nxin] of input samples (see notes) */,
		  int nxout           /* number of x samples output (fast dimension of zout) */, 
		  float dxout         /* x sampling interval output */, 
		  float fxout         /* first x sample output */,
		  int nyout           /* number of y samples output (slow dimension of zout) */, 
		  float dyout         /* y sampling interval output */, 
		  float fyout         /* first y sample output */, 
		  unsigned char *zout /* [nyout][nxout] output samples (see notes) */)
/*< blocky interpolation of a 2-D array of bytes 
  Notes:
  The arrays zin and zout must passed as pointers to the first element of
  a two-dimensional contiguous array of unsigned char values.
  
  Constant extrapolation of zin is used to compute zout for
  output x and y outside the range of input x and y. >*/
/******************************************************************************
Author:  James Gunning, CSIRO Petroleum 1999. Hacked from
intl2b() by Dave Hale, Colorado School of Mines, c. 1989-1991
*****************************************************************************/
{         
    int ixout,iyout,iin,jin;
    float xoff,yoff;
	
    xoff=fxout+0.5*dxin-fxin;
    yoff=fyout+0.5*dyin-fyin;
    for (iyout=0;iyout<nyout;iyout++) {
	jin=(int)((iyout*dyout+yoff)/dyin);
	jin=SF_MIN(nyin-1,SF_MAX(jin,0));						
	for (ixout=0;ixout<nxout;ixout++) {
	    iin=(int)((ixout*dxout+xoff)/dxin);
	    iin=SF_MIN(nxin-1,SF_MAX(iin,0));	
	    zout[nxout*iyout+ixout]=zin[nxin*jin+iin];
	}
    }
}
 
#ifdef TEST
main()
{
    int nxin,nyin,nxout,nyout,ixout,iyout;
    float dxin,fxin,dyin,fyin,fxout,dxout,fyout,dyout;
    unsigned char zin[2][2],zout[4][4];
	
    zin[0][0] = 41;		zin[0][1] = 99;
    zin[1][0] = 99;		zin[1][1] = 99;
    nxin=2;  dxin=1.0;  fxin=0.0;
    nyin=2;  dyin=1.0;  fyin=0.0;
    nxout=4;  dxout=dxin*(nxin-1)/(nxout-1);  fxout=0.0;
    nyout=4;  dyout=dyin*(nyin-1)/(nyout-1);  fyout=0.0;
    intl2b(nxin,dxin,fxin,nyin,dyin,fyin,&zin[0][0],
	   nxout,dxout,fxout,nyout,dyout,fyout,&zout[0][0]);
    for (iyout=0; iyout<nyout; iyout++)
	for (ixout=0; ixout<nxout; ixout++)
	    printf("zout[%d][%d] = %d\n",
		   iyout,ixout,zout[iyout][ixout]);
}
#endif
