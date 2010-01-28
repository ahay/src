/*
 * NUM_FONTS is the total number of fonts.
 * Font 0 should be the good old pen font.
 * The Marker routine assumes that fonts 15 and 16
 * are math and misc fonts.
 * Other than that, I suppose you can reorder them as you like
 * but probably it's a good idea to leave the first 17 fonts just
 * like I've got them here.
 *
 * Additional device-independent fonts can be added quite easily.
 * Fonts 100 and up are reserved for hardware device-dependent fonts,
 * however!
 *
 * The ERRFONT (normally font 0; it is defined in gentext.c) MUST be
 * loaded at compile time. (If your error font itself has an error,
 * woe to you!)
 *
 * - Joe Dellinger
 */
#define NUM_FONTS	19

/*
 * Modifications:
 *
 * Hector Urdaneta (SEP), June 26 1996
 *	header files where change from include "..." to include <...>, 
 *      there path is now defined at compile time
 */

/*
 * Here are the fonts that we want included at compile time into gentext.c.
 * All other fonts will be loaded at runtime on demand.
 */
#include <cyrilc.include>
#include <misc.include>
#include <gothgbt.include>
#include <pen.include>
#include <gothgrt.include>
#include <romanc.include>
#include <gothitt.include>
#include <romand.include>
#include <greekc.include>
#include <romans.include>
#include <greeks.include>
#include <romant.include>
#include <italicc.include>
#include <scriptc.include>
#include <italict.include>
#include <scripts.include>
#include <math.include>
#include <hiragana.include>
#include <katakana.include>


#define	BOTTOM	0
#define	BASE	1
#define	HALF	2
#define	CAP	3
#define	TOP	4
#define	LETTER	5
#define	LINE	6
#define	SPACE	7
#define	START	8
#define	END	9

#define NOT_LOADED 0,0,0,0,0,0,0

typedef struct {
	int load;
	char name[10];
	int *dim;
	int   *saddr;
	int *swidthl;
	int *swidthr;
	int *symbol;
	unsigned int *svec;
	int *lig;
} GLYPH;

/*
 * Here are two examples of fonts loaded at compile time.
 * You have to specify all the entries in the GLYPH structure.
 * The "1" in the first slot specifies that this font is already loaded.
 */
GLYPH font[NUMGENFONT] = {

{1, "pen", pen_dim, pen_addr, pen_widthl, pen_widthr, pen_symbol, pen_vec, pen_lig},


{1, "romans", romans_dim, romans_addr, romans_widthl, romans_widthr, romans_symbol, romans_vec, romans_lig},

/*
 * Here's an example of a font loaded at run time.
 * You put in ZERO for the first entry in the GLYPH structure.
 * This specifies that this font is not loaded.
 * You still put the font name in the second slot as for a loaded font.
 * All the rest of the slots should have zeroes put in them
 * to make the compiler happy.
 * The define NOT_LOADED puts in the correct number of zeroes for you.
 */
{0, "romand", NOT_LOADED},

{1, "romanc", romanc_dim, romanc_addr, romanc_widthl, romanc_widthr, romanc_symbol, romanc_vec, romanc_lig},

{0, "romant", NOT_LOADED},

{0, "italicc", NOT_LOADED},

{0, "italict", NOT_LOADED},

{0, "scripts", NOT_LOADED},

{0, "scriptc", NOT_LOADED},

{0, "greeks", NOT_LOADED},

{0, "greekc", NOT_LOADED},

{0, "cyrilc", NOT_LOADED},

{0, "gothgbt", NOT_LOADED},

{0, "gothgrt", NOT_LOADED},

{0, "gothitt", NOT_LOADED},

{0, "math", NOT_LOADED},

{0, "misc", NOT_LOADED},

{0, "hiragana", NOT_LOADED},

{0, "katakana", NOT_LOADED},

};
