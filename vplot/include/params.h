/*
 * All machine dependencies should go in here!
 */
#define STYLE STANDARD
#ifndef PEN_SPOOL
#define PEN_SPOOL  "/tmp"
#endif

/*
 * The path in the next line usually will end with a slash ('/')!
 * This directory can also be set at runtime via the environmental variable
 * VPLOTFONTDIR. The environmental variable, if set, will override the #define.
 */

#ifndef SYSTEM_FONT_DIRECTORY
#define SYSTEM_FONT_DIRECTORY "/usr/local/lib/vplotfonts/"
#endif

#define DEFAULT_COLOR	 7
#define DEFAULT_LINESTYLE 0
#define DEFAULT_FONT	0
#define DEFAULT_HARDCOPY_FONT	3
#define DEFAULT_PREC	0
#define DEFAULT_HARDCOPY_PREC	2

/*
 * Ratio of fatness to size for symbols
 * (ie, how fat will a 1 inch tall symbol's strokes be, in inches?)
 */
#define SYMBFATRATIO	.05

/*
 * These two get re-allocated larger as needed automatically in dovplot.c,
 * but it's a good idea to start out large enough that that doesn't have
 * to be done very often
 */
#define TXBUFLEN 250	/* Initial max length of a text string */
#define VXBUFLEN 250	/* Initial max number of vertices in a polygon */

/*
 * This is the factor we scale our path and up vectors by before
 * running them through the local text coordinate transformation.
 * (The coordinate transformation, being in terms of device units,
 * gets done in integers. If we don't scale up we get severe roundoff
 * problems for small text sizes at odd angles. We can't make this
 * factor too big, though, or we risk very large text overflowing
 * the maximum possible integer.)
 */
#define TEXTVECSCALE	10.

#define FONTCHECK -1990 /* Magic number that identifies font files */
#define NHATCH	20 	/* Maximum number of lines in a hatch pattern */
/*old
#define MAX_COL 511
#define NPAT	512	
*/

#define MAX_COL 32767	/* Maximum color table number, maximum color number */
#define NPAT	32768	/* Maximum internal pattern number, MAX_COL + 1 */
#define MAXDASH 10	/* Maximum number of dot-dashes */
#define MAX_GUN	255	/* Maximum color gun strength */
#define NUMGENFONT 100  /* Number of fonts reserved for gentext fonts */
#define MAXIN	2048	/* Maximum number of input plot files */
#define MAXFLEN	120	/* Maximum length of a file name */
