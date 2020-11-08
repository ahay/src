/* Convert TIFF image to byte RSF. 

Takes: < file.tiff 
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
#include <tiffio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_file out;
    FILE *tiffin;
    TIFF *tiffout;
    char *tiffname, buf[BUFSIZ];
    unsigned char *grey;
    int n1, i2, n2, nbuf;
    short nc;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_settype(out,SF_UCHAR);

    fclose(sf_tempfile(&tiffname,"w"));
    tiffin = fopen(tiffname,"wb");

    while (1) {
	nbuf = fread(buf,1,BUFSIZ,stdin);
	if (nbuf <= 0) break;
	fwrite(buf,1,nbuf,tiffin);
    }

    fclose(tiffin);
    tiffout = TIFFOpen(tiffname,"rb");

    /* Find the width and height of the image */
    TIFFGetField(tiffout, TIFFTAG_IMAGEWIDTH, &n1);
    TIFFGetField(tiffout, TIFFTAG_IMAGELENGTH, &n2);
    TIFFGetField(tiffout,TIFFTAG_SAMPLESPERPIXEL,&nc);

    if (1==nc) {
	sf_putint(out,"n1",n1);
	sf_putint(out,"n2",n2);
    } else {
	sf_putint(out,"n1",nc);
	sf_putint(out,"n2",n1);
	sf_putint(out,"n3",n2);
    }
  
    grey = sf_ucharalloc(n1*nc);

    for (i2 = 0; i2 < n2; i2++) {
	if (TIFFReadScanline(tiffout, grey, i2, 0) < 0) 
	    sf_error("Trouble reading TIFF file");
	sf_ucharwrite(grey,n1*nc,out);
    }

    TIFFClose(tiffout);
    unlink(tiffname);

    exit(0);
}
