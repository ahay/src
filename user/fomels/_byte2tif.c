/* Convert byte RSF to a TIFF image. */
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
#include <stdio.h>
#include <tiffio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    bool color;
    int n1, i2, n2, nc, nbuf;
    unsigned char *grey=NULL;
    sf_file in=NULL;
    TIFF *tiffout=NULL;
    FILE *tiffin=NULL;
    char *tiffname=NULL, buf[BUFSIZ];

    sf_init(argc,argv);
    in = sf_input("in");

    fclose(sf_tempfile(&tiffname,"w"));
    tiffout = TIFFOpen(tiffname,"wb");

    if (tiffout == NULL) sf_error("can't open file %s\n", tiffname);

    if (SF_UCHAR != sf_gettype(in)) sf_error("Need unsigned char in input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_getbool("color",&color)) color=(bool)(3==n1);

    if (color) {
	nc = n1;
	if (!sf_histint(in,"n2",&n1)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n2)) sf_error("No n3= in input");
    } else {
	nc = 1;
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    }

    grey = sf_ucharalloc (n1*n2*nc);
    sf_ucharread(grey,n1*n2*nc,in);    

    TIFFSetField(tiffout,TIFFTAG_IMAGEWIDTH,n1);
    TIFFSetField(tiffout,TIFFTAG_IMAGELENGTH,n2);
    TIFFSetField(tiffout,TIFFTAG_SAMPLESPERPIXEL,nc);
    TIFFSetField(tiffout,TIFFTAG_BITSPERSAMPLE,8);
    TIFFSetField(tiffout,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
    TIFFSetField(tiffout,TIFFTAG_PHOTOMETRIC, 
		 (3==nc)? PHOTOMETRIC_RGB: PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tiffout,TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiffout, TIFFTAG_ROWSPERSTRIP, 1);

    for (i2 = 0; i2 < n2; i2++) {
	if (TIFFWriteScanline(tiffout, grey+i2*n1*nc, i2, 0) < 0) 
	    sf_error("Trouble writing TIFF file");
    }

    TIFFClose(tiffout);
    tiffin = fopen(tiffname,"rb");

    while (1) {
	nbuf = fread(buf,1,BUFSIZ,tiffin);
	if (nbuf <= 0) break;
	fwrite(buf,1,nbuf,stdout);
    }

    fclose(tiffin);
    unlink(tiffname);

    exit(0);
}
