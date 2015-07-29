/* include files */
/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "../include/extern.h"

/* inline coordinate transforms;  Xpen origin- upper left, pen origin- lower left */
#define	XCORD(x)	x
#define	YCORD(y)	(dev.ymax-y-TBORDER)

/* global definitions */
#define	TBORDER 0
#define	BORDER 0

extern Display *pen_display;
extern int	pen_screen;
extern Window	pen_win;
extern GC	pen_gc;
extern Colormap pen_colormap;

/* pixel map for color manipulation */
extern unsigned long map[256];
/* do we own the whole colormap ? */
extern int own_colormap;
