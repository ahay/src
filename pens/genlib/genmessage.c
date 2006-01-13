/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/genlib/genmessage.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include	<stdio.h>
#include	"../include/mesgcom.h"
#include	"../include/enum.h"

void genmessage (int command, char *string)
/*< Device independent subroutine to handle message operations >*/
{

    switch (command)
    {
    case MESG_HOME:
    case MESG_READY:
	fflush (stderr);
	break;
    case MESG_TEXT:
	fprintf (stderr, "%s", string);
	break;
    case MESG_HIGHLIGHT_ON:
	/* Beep at them to get their attention */
	fprintf (stderr, "\07\07");
	break;
    case MESG_ON:
    case MESG_OFF:
    case MESG_ERASE:
    case MESG_HIGHLIGHT_OFF:
    case MESG_DONE:
    default:
	fflush (stderr);
	break;
    }
}
