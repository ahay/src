#define PS_BLACK		0
#define PS_WHITE		1
#define PSPERIN			72    /* PostScript units per inch */
#define CMPERIN			2.54  /* Centimeters per inch */
#define DEFAULT_BORDER		.25   /* Default unplottable border .25 inch */
#define DEFAULT_PIXELS_PER_INCH	300.
/* 
 * The Apple LaserWriter has a limit of 1500 points for ALL
 * active paths, including the clipping path. Pspen checks only
 * the current drawing path, so we use a lower limit. When this limit
 * is reached, the current path is ended and a new one is begun.
 */
#define PATHLENGTH		1400  

