/*
 * command passed to dev.erase()
 */
#define ERASE_START		0
#define ERASE_MIDDLE		1
#define ERASE_END		2
#define ERASE_BREAK		3

/*
 * Types of erases: (Who would believe that there are 4 kinds?)
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
 * For most screen devices, ERASE_START and ERASE_MIDDLE will just be regular
 * erases, and ERASE_END and ERASE_BREAK will be ignored.
 *
 * Hardcopy devices may use ERASE_MIDDLE and ERASE_END instead.
 */
