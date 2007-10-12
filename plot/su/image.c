/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2001.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 1998.*/
/* All rights reserved.                       */

#include <X11/Xlib.h>
/*^*/

#include <rsf.h>

#include "colormap.h"

XImage *xNewImage (Display *dpy         /* display pointer */, 
		   unsigned long pmin   /* minimum pixel value (corresponding to byte=0) */, 
		   unsigned long pmax   /* maximum pixel value (corresponding to byte=255) */,
		   int width            /* number of bytes in x dimension */, 
		   int height           /* number of bytes in y dimension */, 
		   float blank          /* portion for blanking (0 to 1) */, 
		   unsigned char *bytes /* unsigned bytes to be mapped to an image */)
/*< make a new image of pixels from bytes >*/
/******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 06/08/90
*****************************************************************************/
{
	int scr=DefaultScreen(dpy);
	int i,j,k,line,iline,jline,widthpad;
	float base,scale;
	unsigned long map[256],bkgnd;
	unsigned char *data;
	int byte_perpixel;
	unsigned int depth;
	XImage *xim;
	
	xim=(XImage *) NULL;

	depth=(unsigned int)DefaultDepth(dpy,scr);	
	byte_perpixel=4;
	if(depth<=8) byte_perpixel=1;
	else if(depth<=16) byte_perpixel=2;

/*	else if(depth<=24) byte_perpixel=3;*/



	/* build map for translating bytes to pixels */
	base = ((double) pmin)+0.499;
	scale = ((double) (pmax-pmin))/255.0;
	for (i=0; i<=255; ++i){
		map[i] = base+i*scale;
	}

	/* blanking */
	bkgnd = (unsigned long) WhitePixel(dpy,scr);
	j = SF_MAX(0,SF_MIN(256,(int)(256*blank)));
	for (i = 0; i < j; i++)
		map[255-i] = bkgnd;

	/* allocate memory for image data */
	widthpad = (1+(width-1)/(BitmapPad(dpy)/8))*BitmapPad(dpy)/8;
	data = sf_alloc(widthpad*height,byte_perpixel);

	xim=XCreateImage(	(Display *) dpy,
				(Visual *) DefaultVisual(dpy,scr),
				(unsigned int) DefaultDepth(dpy,scr),
				(int) ZPixmap,
				(int) 0,
				(char *) data,
				(unsigned int) widthpad,
				(unsigned int) height,
/*				(int) BitmapPad(dpy),
				(int) widthpad*byte_perpixel
*/
				8,0);
	byte_perpixel=xim->bits_per_pixel/8;
/*	fprintf(stderr,"\nbyte_perpixel = %d, depth= %d\n", byte_perpixel,depth); */

	/* translate bytes to pixels, padding scanlines as necessary */
	for (line=0; line<height; line++) {
		iline = line*width;
		jline = line*widthpad;
		for (i=iline,j=jline,k=0; k<width; ++i,++j,++k)
		{       if(byte_perpixel==1)
			((unsigned char *)data)[j] =(unsigned char)map[bytes[i]];
			if(byte_perpixel==2)
			  {
			    int edn=xim->byte_order;
			    if(edn==LSBFirst){
			      ((unsigned char *)data)[j*2+0] =(unsigned char)(truecolor_pixel[bytes[i]]);
			      ((unsigned char *)data)[j*2+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>8);
			    }else{
			      ((unsigned char *)data)[j*2+0] =(unsigned char)(truecolor_pixel[bytes[i]]>>24); 
			      ((unsigned char *)data)[j*2+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>16);
			      }

			    /*((unsigned short *)data)[j] =(unsigned short)(truecolor_pixel[bytes[i]]);*/
			  }	
                        if(byte_perpixel==3){
			  int edn=xim->byte_order;
			  if(edn==LSBFirst){
			    ((unsigned char *)data)[j*3+0] =(unsigned char)(truecolor_pixel[bytes[i]]);
			    ((unsigned char *)data)[j*3+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>8);
			    ((unsigned char *)data)[j*3+2] =(unsigned char)(truecolor_pixel[bytes[i]]>>16);
			  }else{
			    ((unsigned char *)data)[j*3+0] =(unsigned char)(truecolor_pixel[bytes[i]]>>24);   
			    ((unsigned char *)data)[j*3+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>16); 
			    ((unsigned char *)data)[j*3+2] =(unsigned char)(truecolor_pixel[bytes[i]]>>8);
			  }
			  
			}	
			if(byte_perpixel==4){
			  int edn=xim->byte_order;
			  if(edn==LSBFirst){
			    ((unsigned char *)data)[j*4+0] =(unsigned char)(truecolor_pixel[bytes[i]]);
			    ((unsigned char *)data)[j*4+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>8);
			    ((unsigned char *)data)[j*4+2] =(unsigned char)(truecolor_pixel[bytes[i]]>>16);
			    ((unsigned char *)data)[j*4+3] =(unsigned char)(truecolor_pixel[bytes[i]]>>24);
			  }else{
			    ((unsigned char *)data)[j*4+0] =(unsigned char)(truecolor_pixel[bytes[i]]>>24);
			    ((unsigned char *)data)[j*4+1] =(unsigned char)(truecolor_pixel[bytes[i]]>>16);   
			    ((unsigned char *)data)[j*4+2] =(unsigned char)(truecolor_pixel[bytes[i]]>>8); 
			    ((unsigned char *)data)[j*4+3] =(unsigned char)(truecolor_pixel[bytes[i]]);
			  }
			    /* ((unsigned long *)data)[j] =(unsigned long)truecolor_pixel[bytes[i]];*/
			}

		}
		for (j=jline+width,k=width; k<widthpad; ++j,++k)
		{

		       if(byte_perpixel==1)
                        ((unsigned char *)data)[j] =((unsigned char *)data)[jline+width-1];
                        if(byte_perpixel==2)
			{
                       /* ((unsigned short *)data)[j] =((unsigned short *)data)[jline+width-1];*/
			 ((unsigned char *)data)[j*2+0] =((unsigned char *)data)[(jline+width-1)*2+0];
                        ((unsigned char *)data)[j*2+1] =((unsigned char *)data)[(jline+width-1)*2+1];
			}
                        if(byte_perpixel==3){
                        ((unsigned char *)data)[j*3+0] =((unsigned char *)data)[(jline+width-1)*3+0];
                        ((unsigned char *)data)[j*3+1] =((unsigned char *)data)[(jline+width-1)*3+1];
                        ((unsigned char *)data)[j*3+2] =((unsigned char *)data)[(jline+width-1)*3+2];
			}
                        if(byte_perpixel==4)
			{
                       ((unsigned char *)data)[j*4+0] =((unsigned char *)data)[(jline+width-1)*4+0];
                        ((unsigned char *)data)[j*4+1] =((unsigned char *)data)[(jline+width-1)*4+1];
                        ((unsigned char *)data)[j*4+2] =((unsigned char *)data)[(jline+width-1)*4+2];
                        ((unsigned char *)data)[j*4+3] =((unsigned char *)data)[(jline+width-1)*4+3];    
			}
		}
	}
	
	/* create and return image structure */
	return xim;


}




