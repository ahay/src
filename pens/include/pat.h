/*
 * area fill pattern stuff
 */

struct pat
{
	int    ydim;
	int    xdim;
	int    ydim_orig;
	int    xdim_orig;
	int   *patbits;
};

extern struct pat pat[];
