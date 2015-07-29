
#define F_COL		(0001)
#define F_FAT		(0002)
#define F_CLIP		(0004)
#define F_FONT		(0010)
#define F_JUST		(0020)
#define F_DASH		(0040)
#define F_COLT		(0100)
#define F_OVLY		(0200)

#if 0
#define VPPEN_NUM_COL	(16384)
#endif
#define VPPEN_NUM_COL	(256)
extern int vpscoltabinfo[VPPEN_NUM_COL][4];
#define ISITSET		0

extern bool vpdumb;
extern int vpstat;
extern int vpalign;
extern int vpfit;
extern float xsize, ysize;
extern bool vpbig;
extern bool vpstyle;
extern const char *vpaligns;
extern int vpcolor, vpfat;
extern int vpsetflag;
extern bool vpblast;
extern int vpbit;
extern int vparray[];
extern int vpasize[];
extern int vpframe, vpframecount;

extern int default_hshift, default_vshift;
extern float default_scale, default_xscale, default_yscale, fatmult_orig;

#include <stdlib.h>
#include <string.h>

