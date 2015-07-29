/*
 * commands to the device's "attributes" routine
 */

#define	SET_COLOR		1
#define SET_COLOR_TABLE 	2
#define SET_WINDOW		3
#define NEW_DASH		4
#define NEW_PAT			5
#define NEW_FONT		6
#define NEW_OVERLAY		7
#define NEW_ALIGN		8
#define NEW_FAT			9
#define BEGIN_GROUP		10
#define END_GROUP		11

#define DOVPLOT_CONT		0	/* MUST be zero */
#define DOVPLOT_EXIT		1

#define UNSPECIFIED_EGROUP	0
#define USER_EGROUP		1
#define ERASE_EGROUP		2
#define IGNORE_EGROUP		3
#define BREAK_EGROUP		4
#define EOF_ERASE_EGROUP	5
#define EOF_IGNORE_EGROUP	6

/*
 * Uses:
 * Dummy arguments are not used, but should be listed to make some finicky
 * compilers happy, even if you don't use any of the commands that require
 * all 4 arguments. By the same token, you should explicitly declare a
 * return value of 0 even for those calls that don't use the return value.
 *
 * Calls will be made at the start of each frame of plotting to
 * set the fatness, alignment mode, font, dash-line style, and raster
 * overlay mode, so you don't need to specifically initialize these.
 *
 *
 * dev.attributes(SET_COLOR, col, dummy1, dummy2, dummy3)
 *
 * Set the current drawing color to col. The device is assumed to
 * start out with the current drawing color 7 = WHITE.
 *
 * dev.attributes(SET_COLOR_TABLE, col, red, green, blue)
 * 
 * Set color number col to have the color (red,green,blue).
 * Red, green, and blue are in the range from 0 to MAX_GUN (255),
 * with 0 being off and MAX_GUN being fully on.
 * Thus, (0,0,0) is black, (128,128,128) is grey, (255,255,255) is white,
 * (255,0,0) is red, (255,255,0) is yellow, etc... The device is assumed
 * to start out with colors 0 through 7 set as in the vplot standard.
 *
 * dev.attributes(SET_WINDOW, xmin, ymin, xmax, ymax)
 *
 * Set current clipping window. (xmin,ymin) is the lower-leftmost
 * displayable point. (xmax,ymax) is the upper-rightmost displayable point.
 * Only necessary to support this if smart_clip=YES, smart_raster=YES, or
 * you don't use a generic routine (genarea) to clip polygons for you.
 * After dev.reset has been called, dovplot assumes that the clipping
 * window is set by the device to the edge of the device's screen.
 * Dev.attributes(SET_WINDOW, ...) will thereafter be called whenever
 * the clipping window is changed from its previous setting. (Thus,
 * if no clipping window is ever set at all, dev.attributes(SET_WINDOW, ...)
 * will never be called.) It is possible for xmin > xmax or ymin > ymax.
 * If this happens EVERYTHING SHOULD BE CLIPPED AWAY.
 *
 * dev.attributes(NEW_DASH, dashon, dummy1, dummy2, dummy3)
 *
 * Warn the device that the current dashed line pattern has just
 * been changed. "dashon" is the new value of the external variable
 * dashon. (See the dev.vector documentation.) If the new current
 * dashed line pattern happens to be just a continuous line,
 * dashon will be 0. In particular this will the be the case at
 * the start of each plot frame.
 *
 * dev.attributes(NEW_PAT, ipat, dummy1, dummy2, dummy3)
 *
 * Warn the device that raster pattern "ipat" has been defined
 * or redefined. The device does not have to allow re-definition
 * of defined patterns (except pattern number 0). (See the dev.area
 * documentation for a description of how patterns are stored.)
 * It is possible for the user to ask for a pattern that never
 * was defined. If this happens, the pattern to use is of dimension
 * "xdim" by "ydim", with all pixels color 0 except the last,
 * which is of the current color. In particular, pattern number 0
 * will always fall into this category (this is a VP_OLDAREA polygon).
 * Patterns numbered > 0 in this category will always have
 * xdim = ydim = 1 (this is just a polygon filled solidly with the
 * current color, which is what dovplot will make if you try to fill
 * with a pattern that you haven't defined).
 * Because all different VP_OLDAREA polygons use just one pattern number,
 * pattern 0, this pattern must be allowed to be redefined. It is OK
 * to use some sort of hardware-dependent grey level instead of
 * the actual pattern if that's the best you can do.
 *
 * dev.attributes(NEW_FONT, txfont, txprec, txovly, dummy1)
 *
 * Warn the device that the text font, precision, or overlay may have
 * been changed. A value of "-1" means no change from the previous value,
 * which may be a device-dependent default. This form of dev.attributes
 * will be called whenever the user uses the "set text font and precision"
 * vplot command (and only then), so you probably should check to see if
 * anything really changed before bothering the device. Even better,
 * only check these variables at the time that you output hardware text,
 * and don't bother with this option in dev.attributes at all.
 *
 * dev.attributes(NEW_OVERLAY, overlay, dummy1, dummy2, dummy3)
 *
 * Warn the device that the raster overlay mode may have been changed. Again,
 * most devices will ignore this dev.attributes call, and check the
 * external integer "overlay" only when they actually need to know what
 * it is.
 *
 * dev.attributes(NEW_ALIGN, txalign.hor, txalign.ver, dummy1, dummy2)
 *
 * Warn the device that the text alignment mode has been changed. Again,
 * most devices will ignore this call and just check the alignment
 * when hardware text is drawn.
 *
 * dev.attributes(NEW_FAT, fat, dummy1, dummy2, dummy3)
 *
 * Warn the device that the fatness has been changed. Again, most devices
 * should ignore this call and just check the fatness before processing
 * any of the 3 things in which it is used: vectors, text, markers.
 * Check the dev.vector documentation for caveats similar to those for
 * the NEW_PAT case.
 *
 * dev.attributes(BEGIN_GROUP, group_number, pos, dummy1, dummy2)
 * dev.attributes(END_GROUP, group_number, type, dummy1, dummy2)
 *
 * extern char group_name[MAXFLEN+1];
 *
 * return iflag
 *
 * Everything between these calls is grouped as one "object".
 * Objects may be nested. Objects may not extend across erases,
 * breaks, or files. The name of the most recently begun object is
 * stored in the external variable group_name. Group_name is defined
 * in extern.h. "group_number" gives the heirarchical order of the
 * object. Each plot frame is itself an object of group number 0.
 * "pos" gives the offset from the start of the file of the start of
 * this object. If "pos" is negative, then the input file is one
 * for which seeks don't work. Dev.attributes(END_GROUP,...) returns
 * "iflag" to tell dovplot what it should do next.
 * If "iflag" is DOVPLOT_CONT (=0), then dovplot should continue exactly
 * as normal. If "iflag" is DOVPLOT_EXIT, then dovplot will exit upon the
 * return of this call. Note that K+R says that failing to declare
 * a return value returns garbage, so you MUST declare a return value
 * for dev.attributes(END_GROUP,...). To help make sure you do this,
 * dovplot will issue a warning message if an unknown return value is
 * given.
 *
 * The following esoterica is useful for devices that use dovplot as a
 * "frame plotter", with a frame beginning with a call to dovplot and ending
 * with dovplot's being told to exit after calling
 * dev.attributes(END_GROUP,...).
 * "type" is passed from dovplot to dev.attributes(END_GROUP,...)
 * to inform the device what caused this particular END_GROUP to occur.
 * UNSPECIFIED_EGROUP means that you've been called by an old version of
 * dovplot (or some other subroutine) that didn't know about this argument
 * and just passed 0 intended as a dummy placeholder. In such a case complain,
 * but continue in some reasonable default way.
 * USER_EGROUP means it is an end group that closes a user-defined group
 * (group level 1 or greater).
 * The following classes do NOT correspond to user-defined groups, but are
 * generated within dovplot itself (level-0 groups).
 * ERASE_EGROUP means dovplot just hit an erase (or a break that it has been
 * told to treat as an erase via break=erase on the command line) within the
 * input plotfile, and so is ending the current internally generated level-0
 * group that it started either at the beginning of the current plotfile
 * or after the last erase or break.
 * BREAK_EGROUP is similar to ERASE_EGROUP except dovplot just hit a break
 * not an erase (and breaks are not being ignored or treated as erases).
 * Note that after a break a new frame begins WITHOUT ERASING THE STUFF LEFT
 * ON THE SCREEN FROM THE PREVIOUS FRAME.
 * IGNORE_EGROUP means dovplot just hit an erase or break, but dovplot
 * has been told to "ignore" such things within plotfiles (erase=n on the
 * command line). In such cases dovplot makes all the calls to do prompting,
 * etc, but never gets around to calling dev.erase at all. (Unlike the case
 * when there is an unignored and unchanged break in the input plotfile, when
 * dev.erase IS called with "ERASE_BREAK".)
 * Probably you want to treat IGNORE_EGROUP just like BREAK_EGROUP.
 * EOF_ERASE_EGROUP means dovplot got to the end of the current plotfile,
 * and the command-line erase setting is such that there is an automatic
 * erase at the beginning of processing of each plotfile (erase=yes). You
 * probably want to treat it just like a regular endgroup generated by
 * an erase within a file.
 * EOF_IGNORE_EGROUP means dovplot got to the end of the current plotfile,
 * but the command-line erase setting is such that there will be NO dovplot-
 * generated erase at the beginning of processing of the next plotfile
 * (could be erase=literal or erase=no). Treat this end-of-file like an
 * automatically generated break instead of an automatically generated erase.
 * (Note it is possible that erase=literal and the next plotfile happens to
 * begin with an erase! You may need to keep track of whether dev.erase gets
 * called and told to erase the screen before the level-0 group beginning the
 * next plotfile starts up. Or you may just have to punt, which isn't SO bad
 * because the erase=literal option is only rarely used.)
 *
 */
