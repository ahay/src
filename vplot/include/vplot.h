#ifndef VPLOT_H
#define VPLOT_H
#include<stdio.h>

/*
 * Weird backwards-compatible units
 */
#define RPERIN 		600.	/* vplot units per inch */
#define HATCHPERIN	100.	/* Hatch units per inch */
#define TXPERIN 	33.	/* Text units per inch */
#define FATPERIN	200.	/* Fatness units per inch */
/*
 * Height in inches of "standard" device, standard style
 */
#define STANDARD_HEIGHT 10.24
/*
 * Height in inches of "standard" device, rotated style
 */
#define ROTATED_HEIGHT	 7.5
/*
 * Aspect ratio of the default window (height/width)
 */
#define SCREEN_RATIO 0.75
#define VP_MAX 54.6		/* absolute maximum x or y coordinate in inches */

/*
 * text alignment enumerations
 */
/* horizontal */
#define	TH_NORMAL	0
#define	TH_LEFT		1
#define	TH_CENTER	2
#define	TH_RIGHT	3
#define TH_SYMBOL	4

/* vertical */
#define	TV_NORMAL	0
#define TV_BOTTOM	1
#define TV_BASE		2
#define TV_HALF		3
#define TV_CAP		4
#define TV_TOP		5
#define TV_SYMBOL	6

struct txalign {
	int hor;
	int ver;
};

/*
 * text precision enumerations
 */
#define STRING	0
#define CHAR	1
#define STROKE	2
/* leave it what it already was */
#define NO_CHANGE -1

/*
 * text overlay enumerations
 */
#define OVLY_NORMAL	0
#define OVLY_BOX	1
#define OVLY_SHADE	2
#define OVLY_SHADE_BOX	3

/*
 * colors
 */
#define BLACK    0
#define BLUE     1
#define RED      2
#define PURPLE   3
#define GREEN    4
#define CYAN     5
#define YELLOW   6
#define WHITE    7

/*
 * Coordinate Origin
 */
#define STANDARD	0
#define ROTATED		1
#define ABSOLUTE 	3

/*
 * Fonts
 */

#define PEN		0
#define ROMANS		1
#define ROMAND		2
#define ROMANC		3
#define ROMANT		4
#define ITALICC		5
#define ITALICT		6
#define SCRIPTS		7
#define SCRIPTC		8
#define GREEKS		9
#define GREEKC		10
#define CYRILC		11
#define GOTHGBT		12
#define GOTHGRT		13
#define GOTHITT		14
#define MATH		15
#define MISC		16

/*
 * vplot metafile op-codes
 */

#define VP_SETSTYLE		'S'

#define VP_MOVE			'm'
#define VP_DRAW			'd'
#define VP_PLINE	    	'L'
#define VP_PMARK	   	'M'
#define VP_TEXT			'T'
#define VP_GTEXT		'G'
#define VP_AREA			'A'
#define VP_OLDAREA		'a'
#define VP_BYTE_RASTER		'R'
#define VP_BIT_RASTER		'r'
#define VP_SHORT_RASTER		'B'

#define VP_ERASE		'e'
#define VP_BREAK		'b'
#define VP_PURGE		'p'
#define VP_NOOP			'n'

#define VP_ORIGIN		'o'
#define VP_WINDOW		'w'

#define VP_FAT			'f'
#define VP_SETDASH		's'
#define VP_COLOR		'c'
#define VP_SET_COLOR_TABLE	'C'
#define VP_TXALIGN		'J'
#define VP_TXFONTPREC		'F'
#define VP_PATLOAD		'l'
#define VP_OVERLAY		'v'

#define VP_MESSAGE		'z'
#define VP_BEGIN_GROUP		'['
#define VP_END_GROUP		']'

/* Hopefully now dead primitives */
#define VP_OLDTEXT		't'


extern short geth(register FILE*);
extern int name_to_coltab(char *colname, int nocol, float *red, float *green, float *blue);
extern short puth (register int w, register FILE *iop);
extern int vp_area(float *xp, float *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_uarea(float *xp, float *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_arrow (float x0,float  y0,float  x,float  y,float r);
extern int vp_uarrow (float x0,float  y0,float  x,float  y,float r);
extern int vp_bgroup(char *string);
extern int vp_break();
extern int vp_clip(float, float, float, float);
extern int vp_uclip(float, float, float, float);
extern int vp_color(int col);
extern int vp_coltab(int,float,float,float);
extern int vp_dash(float,float,float,float);
extern int vp_draw(float,float);
extern int vp_udraw(float,float);
extern int vp_egroup();
extern int vp_endplt();
extern int vp_erase();
extern int vp_fat(int FATNESS);
extern int vp_file(char *filename);
extern int vp_filep(FILE *filepntr);
extern int vp_fill(float*,float*,int);
extern int vp_ufill(float*,float*,int);
extern int vp_gtext(float X,float Y,float  XPATH,float  YPATH,float  XUP,float  YUP,char *string);
extern int vp_hatchload(int ANGLE,int  NUMHATCH, int IHATCH, int *hatcharray);
extern int vp_message(char*); 
extern int vp_move(float,float);
extern int vp_umove(float,float);
extern int vp_orig(float,float);
extern int vp_uorig(float,float);
extern int vp_patload(int PPI,int  NX,int  NY,int  IPAT,int  *colarray);
extern int vp_pendn(float,float); 
extern int vp_upendn(float,float); 
extern int vp_penup();
extern int vp_pline(float *XP, float *yp,int LP);
extern int vp_upline(float *XP, float *yp,int LP);
extern int vp_plot(float,float,int);
extern int vp_pmark(int NPTS,int  MTYPE, int MSIZE, float *xp,float *yp);
extern int vp_upmark(int NPTS,int  MTYPE, int MSIZE, float *xp,float *yp);
extern int vp_purge();
extern int vp_rascol16tab(int nreserve,char *colname);
extern int vp_rascoltab(int nreserve,char *colname);
extern int vp_rastershort (unsigned short *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, float XLL, float YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_raster (unsigned char *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, float XLL, float YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_scale(float,float);
extern int vp_setdash(float *dashp,float *gapp, int LP);
extern int vp_stretch(float XMIN,float YMIN,float  XMAX,float YMAX);
extern int vp_style(int);
extern int vp_text(float X,float  Y,int  SIZE,int  ORIENT, char *string);
extern int vp_tfont(int,int,int);
extern int vp_tjust(int,int);
extern int vp_uclip(float,float,float,float);
extern int vp_gtext(float X,float Y, float XPATH,float YPATH,float XUP,float YUP, char *string);
extern int vp_ugtext(float X,float Y, float XPATH,float YPATH,float XUP,float YUP, char *string);
extern int vp_uplot(float X, float Y,int  DOWN);
extern int vp_uraster (unsigned char *array, int BLAST,int  BIT,int  OFFSET, int  XPIX,int  YPIX,float  XLL,float  YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_utext(float X,float Y,int SIZE,int ORIENT, char *string);
extern int vp_where(float*,float*);
extern void p_pout (float xp,float  yp,int  down, FILE *plt);







#endif
