/*
 * commands to the device's "interact" routine
 */

#define	INT_USER_PAUSE		0
#define	INT_GET_STRING		1
#define	INT_PAUSE		2
#define	INT_F_PAUSE		3

/*
 * These are defined in attrcom.h but also used by dev.interact
 * #define DOVPLOT_CONT            0	(this has to be 0)
 * #define DOVPLOT_EXIT            1
 */

/*
 * Termout is conncted to read input from "/dev/tty".
 * You don't have to use it if that isn't appropriate.
 *
 * int dev.interact(INT_GET_STRING, termout, instring)
 * FILE *termout;
 * char *instring;
 *
 * Wait for the user to enter a newline-terminated string.
 * Put the user's string into instring and return. (This option
 * is currently not used, but is provided to allow interactive
 * extensions in the future without changing the device-dependent code
 * interface.)
 *
 *
 * int dev.interact(INT_PAUSE, termout, instring)
 * FILE *termout;
 * char *instring;
 *
 * Wait for the user to indicate his desire for plotting to continue,
 * in whatever fashion is appropriate for that device, then return.
 * For most devices, this will be by hitting "return" on the keyboard,
 * and so this case can simply be handled the same way as the previous one.
 * Any text stuffed into "instring" will be ignored, however.
 *
 * Note dovplot checks the return code returned by dev.interact in the
 * INT_PAUSE case. If it is DOVPLOT_EXIT, then dovplot will return right away
 * instead of continuing to process the next frame in the plotfile
 * (if there is a next frame). Normally dev.interact should return
 * DOVPLOT_CONT (0).
 *
 *
 * int dev.interact(INT_F_PAUSE, termout, instring)
 * FILE *termout;
 * char *instring;
 *
 * Just like the INT_PAUSE case, except that this is a pause generated
 * by "endpause=y", occuring after all plotting is over with.
 *
 * int dev.interact(INT_USER_PAUSE, termout, instring)
 * FILE *termout;
 * char *instring;
 *
 * Never used from dovplot. Provided so that in case the device wants to
 * make its own special case for internal use only, a number is already
 * reserved.
 */
