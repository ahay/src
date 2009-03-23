extern unsigned char *image;	/* raster image */
extern int rascolor;		/* current color */
extern int esize;
#define	Image(IX,IY)		image[(IX)+dev.xmax*(dev.ymax-1-(IY))]
#define	Image3(IX,IY,RGB)	image[(RGB)+(IX)*3+dev.xmax*3*(dev.ymax-1-(IY))]
#define	Min(IX,IY)		((IX) < (IY) ? (IX) : (IY))
#define	Max(IX,IY)		((IX) > (IY) ? (IX) : (IY))
#define NCOLOR 256		/* number of colors */
