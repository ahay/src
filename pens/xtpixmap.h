/*
 * Copyright 1991 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/xtlib/xtpixmap.h
 *
 * Steve Cole (SEP), August 4 1991
 *	Wrote xtpen.  
 *	Inserted this sample edit history entry.
 */

extern Pixmap pen_pixmap;
extern int pen_width,pen_height;

extern int have_pixmap;

extern void clear_pixmap();
extern void create_pixmap();
extern void remove_pixmap();

#if defined(__STDC__) || defined (__stdc__)
extern Pixmap MyCreatePixmap( Display *, Drawable,unsigned int,
		unsigned int,unsigned int );
#else
extern Pixmap MyCreatePixmap();
#endif
