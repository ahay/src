/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*
 * bheader.h - include file for binary header sizes
 * THIS HEADER FILE IS FIXED FOR ALL MACHINES BY THE SEG_Y STANDARD 
 */

#ifndef BHEADER_H
#define BHEADER_H

#define	EBCBYTES	3200	/* Bytes in the card image EBCDIC block */
#define	BNYBYTES	400	/* Bytes in the binary coded block	*/
#define	SEGY_HDRBYTES	240	/* Bytes in the tape trace header	*/
#define	SEGY_NKEYS	71	/* Number of mandated header fields	*/
#define BHED_NKEYS	27	/* Number of mandated binary fields	*/

#endif
