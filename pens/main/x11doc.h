
const char *documentation[] = {
"NAME",
#ifdef SEP
"	X11pen - SEPLIB vplot filter for X-Windows",
#else
"	x11pen - vplot filter for X-windows",
#endif
"",
"SYNOPSIS",
#ifdef SEP
"	X11pen [options] in=vplot-input file OR header file OR stdin",
#else
"	x11pen [options] [inputfiles]",
#endif
"",
"OPTIONS",
"    X=[host,remote-host] (X-server)  |  stay=n (y keeps image after exit)",
"    width=.77(screen)  | height=.56(screen) | align=tr(top,mid,bot;lef,cen,rig)",
#include "../include/gendoc.h"
"",
"SEE ALSO",
"	man pen"
};
int doclength = {sizeof documentation/sizeof documentation[0]};

