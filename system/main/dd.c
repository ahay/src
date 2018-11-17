/* Convert between different formats. */
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
#include <rsf.h>

static void ddbreak (sf_datatype itype, sf_datatype otype);
static float ibm2float (const char* num);

int main(int argc, char *argv[])
{
    off_t size;
    size_t nin, nout, bufsiz, ein, eout;
    int line, n1, i, j, *ibuf, strip;
    sf_file in, out;
    char *form, *type, *format, *bufin, *bufout;
    sf_datatype itype, otype;
    float *fbuf;
    short *sbuf;
    off_t *lbuf;
    sf_complex *cbuf;
    unsigned char *ubuf;
    bool ibm = false;
    bool trunc = false;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output ("out");
    size = sf_filesize(in);
    if (!sf_getbool("trunc", &trunc)) trunc=false;
    /* Truncate or round to nearest when converting from float to int/short */

    form = sf_getstring("form");
    /* ascii, native, xdr */
    type = sf_getstring("type");
    /* int, float, complex, short, long */
    if (NULL == form && NULL == type) sf_error("Specify form= or type=");
	
    if (NULL != form) {
	switch (form[0]) {
	    case 'a':
		sf_setform(out,SF_ASCII);
		if (!sf_getint("line",&line)) line=8;
		/* Number of numbers per line (for conversion to ASCII) */
		format = sf_getstring("format");
		/* Element format (for conversion to ASCII) */
		if (!sf_getint("strip",&strip)) strip=0;
		/* If strip characters from format at the end of the line */
		sf_setaformat(format,line,strip);
		break;
	    case 'n':
		sf_setform(out,SF_NATIVE);
		break;
	    case 'x':
		sf_setform(out,SF_XDR);
		break;
	    default:
		sf_error("Unsupported form=\"%s\"",form);
		break;
	}
    }

    itype = sf_gettype (in);
    otype = itype;
    if (NULL != type) {
	switch (type[0]) {
	    case 'i':
		otype = SF_INT;
		break;
	    case 'f':
		otype = SF_FLOAT;
		break;
	    case 'c':
		otype = SF_COMPLEX;
		break;
            case 's':
                otype = SF_SHORT;
		break;
	    case 'l':
		otype = SF_LONG;
		break;
	    case 'd':
		otype = SF_DOUBLE;
		break;
	    default:
		sf_error("Unsupported type=\"%s\"",type);
		break;
	}
    } 
    sf_settype(out,otype);

    ein = sf_esize(in);
    eout = sf_esize(out);

    /* optimize buffer size */
    bufsiz = sf_bufsiz(in);
    bufin  = sf_charalloc(bufsiz);
    bufout = sf_charalloc(bufsiz);
    bufsiz /= ein;

    if (SF_UCHAR != itype && SF_CHAR != itype && 
	SF_SHORT != itype && SF_SHORT != otype &&
	ein != eout && sf_histint(in,"n1",&n1)) {
	n1 = (n1*ein)/eout;
	sf_putint(out,"n1",n1);
    }

    while (size > 0) {
	nin = (bufsiz < size)? bufsiz:size;
	nout = nin*ein/eout;
	switch (itype) {
	    case SF_LONG:
		nin = bufsiz*ein/eout;
		if (nin > size) nin=size;

		lbuf = (off_t*) bufin;
		sf_longread(lbuf,nin,in);
		switch (otype) {
                    case SF_SHORT:
			sbuf = (short*) bufout;
			for (i=0; i < (int) nin; i++) {
			    sbuf[i] = lbuf[i]; 
			}
                        sf_shortwrite(sbuf,nin,out);
			break;
		    case SF_INT:
                        ibuf = (int*) bufout;
                        for (i=0; i < (int) nin; i++) {
			    ibuf[i] = lbuf[i]; 
			}
                        sf_intwrite(ibuf, nin, out);
			break;
		    case SF_FLOAT:
			fbuf = (float*) bufout;
			for (i=0; i < (int) nin; i++) {
			    fbuf[i] = lbuf[i]; 
			}
			sf_floatwrite(fbuf,nin,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
            case SF_SHORT:
		nin = bufsiz*ein/eout;
		if (nin > size) nin=size;

		sbuf = (short*) bufin;
		sf_shortread(sbuf,nin,in);
		switch (otype) {
                    case SF_SHORT:
                        sf_shortwrite(sbuf,nin,out);
			break;
		    case SF_INT:
                        ibuf = (int*) bufout;
                        for (i=0; i < (int) nin; i++) {
			    ibuf[i] = sbuf[i]; 
			}
                        sf_intwrite(ibuf, nin, out);
			break;
		    case SF_FLOAT:
			fbuf = (float*) bufout;
			for (i=0; i < (int) nin; i++) {
			    fbuf[i] = sbuf[i]; 
			}
			sf_floatwrite(fbuf,nin,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_INT:
		ibuf = (int*) bufin;
		sf_intread(ibuf,nin,in);
		switch (otype) {
		    case SF_INT:
			sf_intwrite(ibuf,nout,out);
			break;
                    case SF_SHORT:
			nout = nin;
                        sbuf = (short*) bufout;
                        for (i=0; i < (int) nin; i++) {
			    sbuf[i] = (short) ibuf[i]; 
			}
                        sf_shortwrite(sbuf, nin, out);
			break;
		    case SF_FLOAT:
			fbuf = (float*) bufout;
                        if (!sf_getbool("ibm",&ibm)) ibm=false;
                        /* Special case - assume integers actually represent IBM floats */
                        if (ibm) {
                            /* This is a very special case - the input values
                               are not integers but IBM floats actually. Since
                               there is no format for storing float headers in Madagascar,
                               they come in as integers. We catch them here and
                               do a proper conversion from IBM floats to IEEE floats */
			    for (i=j=0; i < (int) nin && j < nout; i++, j++) {
			        fbuf[j] = ibm2float ((const char*)&ibuf[i]);
			    }
                        } else {
                            /* Regular conversion from integer to float */
			    for (i=j=0; i < (int) nin && j < nout; i++, j++) {
			        fbuf[j] = ibuf[i]; 
			    }
                        }
			sf_floatwrite(fbuf,nout,out);
			break;
		    case SF_COMPLEX:
			cbuf = (sf_complex*) bufout;
			for (i=j=0; i < (int) nin && j < nout; i+=2, j++) {
			    cbuf[j] = sf_cmplx(ibuf[i],ibuf[i+1]); 
			}
			sf_complexwrite(cbuf,nout,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_FLOAT:
		fbuf = (float*) bufin;
		sf_floatread(fbuf,nin,in);
		switch (otype) {
		    case SF_FLOAT:
			sf_floatwrite(fbuf,nout,out);
			break;
		    case SF_INT:
			ibuf = (int*) bufout;
			/* Avoiding a conditional inside a loop. More code */
			if (trunc) {
			    for (i=j=0; i < nout && j < (int) nin; i++, j++) {
			        ibuf[i] = fbuf[j];
			    }
			} else {
			    for (i=j=0; i < nout && j < (int) nin; i++, j++) {
			        ibuf[i] = roundf(fbuf[j]);
			    }
			}
			sf_intwrite(ibuf,nout,out);
			break;
                    case SF_SHORT:
                        sbuf = (short*) bufout;
			/* Strange why the structure of the loop is different
			   for short than for int. Different esize? Should
			   investigate later */
			if (trunc) {
                            for (i=0; i < (int) nin; i++) {
			        sbuf[i] = (short) fbuf[i];
			    }
			} else {
			    for (i=0; i < (int) nin; i++) {
			        sbuf[i] = roundf(fbuf[i]);
			    }
			}
                        sf_shortwrite(sbuf, nin, out);
			break;
		    case SF_COMPLEX:
			cbuf = (sf_complex*) bufout;
			for (i=j=0; i < nout && j < (int) nin; i++, j+=2) {
			    cbuf[i] = sf_cmplx(fbuf[j],fbuf[j+1]); 
			}
			sf_complexwrite(cbuf,nout,out);
			break;		
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_COMPLEX:
		cbuf = (sf_complex*) bufin;
		sf_complexread(cbuf,nin,in);
		switch (otype) {
		    case SF_COMPLEX:
			sf_complexwrite(cbuf,nout,out);
			break;
		    case SF_FLOAT:
			fbuf = (float*) bufout;
			for (i=j=0; i < nout && j < (int) nin; i+=2, j++) {
			    fbuf[i]   = crealf(cbuf[j]); 
			    fbuf[i+1] = cimagf(cbuf[j]);
			}
			sf_floatwrite(fbuf,nout,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_UCHAR:
		nin = bufsiz*ein/eout;
		if (nin > size) nin=size;

		ubuf = (unsigned char*) bufin;
		sf_charread(bufin,nin,in);
		switch (otype) {
		    case SF_INT:
			ibuf = (int*) bufout;
			for (i=0; i < (int) nin; i++) {
			    ibuf[i] = (int) ubuf[i];
			}
			sf_intwrite(ibuf,nin,out);
			break;
		    case SF_FLOAT:
			fbuf = (float *) bufout;
			for (i=0; i < (int) nin; i++) {
			    fbuf[i] = (float) ubuf[i];
			}
			sf_floatwrite(fbuf,nin,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    case SF_CHAR:
		nin = bufsiz*ein/eout;
		if (nin > size) nin=size;
		sf_charread(bufin,nin,in);
		switch (otype) {
		    case SF_INT:
			ibuf = (int*) bufout;
			for (i=0; i < (int) nin; i++) {
			    ibuf[i] = (int) bufin[i];
			}
			sf_intwrite(ibuf,nin,out);
			break;
		    case SF_FLOAT:
			fbuf = (float *) bufout;
			for (i=0; i < (int) nin; i++) {
			    fbuf[i] = (float) bufin[i];
			}
			sf_floatwrite(fbuf,nin,out);
			break;
		    default:
			ddbreak (itype,otype);
			break;
		}
		break;
	    default:
		ddbreak (itype,otype);
		break;
	}
	size -= nin;
    }
    

    exit(0);
}

static void ddbreak (sf_datatype itype, sf_datatype otype)
{
    const char* types[]={"uchar","char","int","float",
			 "complex","short","double"};

    sf_error("Conversion from %s to %s"
	     " is unsupported",types[itype],types[otype]);
}

static float ibm2float (const char* num)
/* floating point conversion from IBM format */
{
    unsigned int x, s, f;
    const unsigned int fMAXIEEE = 0x7F7FFFFF;
    int e;         
    float y;
                                                                     
    x = *((unsigned int*)num);
    
    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) return 0.0; 

    /* fetch the sign, exponent (removing excess 64), and fraction */   
    s =   x & 0x80000000;                                               
    e = ((x & 0x7f000000) >> 24) - 64;                                   
    f =   x & 0x00ffffff;                                                
                                                                    
    /* convert scale factor from base-16 to base-2 */        
    if (e >= 0) {
	e <<= 2;  
    } else { 
	e = -((-e) << 2); 
    }
                                                                        
    /* convert exponent for 24 bit fraction to 23 bit fraction */           
    e -= 1;                                                               
                                                                            
    /* normalize the fraction */                                            
    if (0 != f) {
	while ((f & 0x00800000) == 0) {         
	    f <<= 1;
	    e -= 1;
	}       
    }                                                               
	
    /* drop the '1' preceeding the binary point */                       
    f &= 0x007fffff;                                                         
    
    /* convert exponent to excess 127 and store the number */
    if ((e += 127) >= 255) {
	s |= fMAXIEEE;
    } else if (e > 0) {
	s |= (e << 23) | f; 	    
    }    

    memcpy (&y,&s,4);
    return y;
}
