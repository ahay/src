/*
 * command passed to dev.erase()
 */
enum {
    ERASE_START, 
    ERASE_MIDDLE, 
    ERASE_END,
    ERASE_BREAK,
    ERASE_BACKGROUND    
    };

/*
 * Types of erases: (Who would believe that there are 5 kinds?)
 *
 * An erase really has two parts; you END one frame of a plot and BEGIN
 * another. For some devices you have to be careful to separate these
 * two different parts.
 *
 * ERASE_START is an erase which occurs before anything has been drawn.
 *
 * ERASE_MIDDLE is the typical case. You are both erasing something
 * 	that is already there, and also preparing for another plot that
 *	is coming.
 *
 * ERASE_END is a terminal erase, with nothing more to be drawn afterwards.
 *
 * ERASE_BREAK is like ERASE_MIDDLE, except that the picture shouldn't
 *	actually be cleared. The purpose of the Break command in Vplot
 *	is to allow saving of snapshots of a plot in progress. (An example
 *	of the proper use of ERASE_BREAK is in the filter Raspen.)
 *
 * ERASE_BACKGROUND is for devices that can't easily change the background
 *	color, for example, a paper plotter. If the vplot file redefines
 *	color 0 to be red, do you really want to color the entire sheet
 *	red? If the user uses the "background" command, it means they do.
 *      Most devices will leave the external variable smart_background
 *      set to its default value of false, in which case dovplot will
 *	honor the vplot background command by drawing a polygon of color 0
 *	that fills the plottable area. There are two cases where you
 *	would set smart_background to true. First, if your device
 *	naturally handles changes to the background color (for example,
 *	a graphics device with a 256-color color table, where changing
 *	color 0 instantly changes what has already been plotted). In
 *	that case the device may safely ignore this command.
 *	Second, if you have a "smart" device that has a built-in
 *	command you can use to change the background color.
 *
 * For most screen devices, ERASE_START and ERASE_MIDDLE will just be regular
 * erases, and ERASE_END and ERASE_BACKGROUND will be ignored.
 *
 * Hardcopy devices will generally use ERASE_MIDDLE and ERASE_END and ignore
 * ERASE_START and ERASE_BREAK.
 */
