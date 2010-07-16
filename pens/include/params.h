/*
 * All machine dependencies should go in here!
 */

/*
 * The path in the next line usually will end with a slash ('/')!
 * This directory can also be set at runtime via the environmental variable
 * VPLOTFONTDIR. The environmental variable, if set, will override the #define.
 */

#ifndef SYSTEM_FONT_DIRECTORY
#define SYSTEM_FONT_DIRECTORY "/usr/local/lib/vplotfonts/"
#endif

#define BACKGROUND_COLOR 0
#define DEFAULT_COLOR	 7
#define DEFAULT_LINESTYLE 0
#define DEFAULT_FONT	0
#define DEFAULT_HARDCOPY_FONT	3
#define DEFAULT_SANSSERIF_FONT	2
#define DEFAULT_PREC	0
#define DEFAULT_HARDCOPY_PREC	2

#define GREEK_SERIF_FONT	10
#define GREEK_SANSSERIF_FONT	9

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

#define FONTCHECK -1990 /* Magic number that identifies font files */
#define NHATCH	20 	/* Maximum number of lines in a hatch pattern */
/*old
#define MAX_COL 511
#define NPAT	512	
*/

#define MAX_COL 32767	/* Maximum color table number, maximum color number */
#define NPAT	32768	/* Maximum internal pattern number, MAX_COL + 1 */
#define MAXDASH 10	/* Maximum number of dot-dashes */
#define NUMGENFONT 100  /* Number of fonts reserved for gentext fonts */
#define MAXIN	2048	/* Maximum number of input plot files */
#define MAXFLEN	120	/* Maximum length of a file name */
