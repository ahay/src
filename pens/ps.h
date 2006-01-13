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
extern int tex;
extern int file_created;
extern int ncopies_document;
extern char label[];
extern char scratch_file[];
extern int hold;
extern int lost;
extern int ps_dash_pattern_set;
extern int ps_dash_pattern_exists;
extern int ps_last_fat;
extern float psscale;
extern int ps_color;
extern int dumb_fat;
extern int force_color;
extern float ps_xlength, ps_ylength;
extern int ps_set_papersize;
extern float ps_ypapersize;
extern char psprintertype[];
