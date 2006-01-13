/*
 * indent destroys this, so I stuffed it down in an include file
 * out of harm's way
 */

#define READ_RASTER(BYTELINE,BYTE2LINE,BITLINE)	\
{\
    if (c == VP_BYTE_RASTER)\
    {\
	for (j = 0; j < num_byte; j++)\
	{\
	    BYTELINE;\
	}\
    }\
    else if (c == VP_SHORT_RASTER)\
    {\
	for (j = 0; j < num_byte; j++)\
	{\
	    BYTE2LINE;\
	}\
    }\
    else\
    {\
	for (jj = 0; jj < num_byte; jj += 8)\
	{\
	    ibyte = (unsigned char) fgetc (pltin);\
	    for (j = 0; j < 8; j++)\
	    {\
		BITLINE;\
	    }\
	}\
    }\
}
