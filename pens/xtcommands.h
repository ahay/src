

/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my eXample below.
 */

/*
 *
 *  source file:   ./filters/xtlib/xtbuttons.c
 *
 * Steve Cole (SEP), February 18 1992
 *      Inserted this sample edit history entry.
 *
 */


extern void actionNext();
extern void actionPrev();
extern void actionQuit();
extern void actionRun();
extern void actionStop();
extern void actionRestart();
extern void actionStretch();
extern void actionSlower();
extern void actionFaster();
extern void actionNumber();
extern void actionNumReset();
extern void actionGoto();
extern void actionCoord();
extern void actionRunMode();

extern void PenRepaint();

/* 	 <ColormapNotify>:	PenRepaint() \n\
         <ConfigureNotify>:	PenRepaint() \n\ */

/* default translation table for pen_picture widget */
static char trans[] =
	"<Expose>:		PenRepaint() \n\
         <ConfigureNotify>:	PenRepaint() \n\
         <Btn1Down>:            xt_print_coord() \n\
         None<KeyPress>n:       xt_stop() xt_reset_number() xt_next() \n\
         None<KeyPress>m:       xt_stop() xt_reset_number() xt_prev() \n\
         None<KeyPress>r:       xt_run()  \n\
         None<KeyPress>q:       xt_quit()  \n\
         None<KeyPress>.:       xt_stop()  \n\
         None<KeyPress>f:       xt_faster()  \n\
         None<KeyPress>s:       xt_slower()  \n\
         None<KeyPress>t:       xt_stretchy()  \n\
	 None<KeyPress>Escape: 	xt_reset_number()\n\
         None<KeyPress>0: 	xt_number(0)\n\
         None<KeyPress>1: 	xt_number(1)\n\
         None<KeyPress>2: 	xt_number(2)\n\
         None<KeyPress>3: 	xt_number(3)\n\
         None<KeyPress>4: 	xt_number(4)\n\
         None<KeyPress>5: 	xt_number(5)\n\
         None<KeyPress>6: 	xt_number(6)\n\
         None<KeyPress>7: 	xt_number(7)\n\
         None<KeyPress>8: 	xt_number(8)\n\
         None<KeyPress>9: 	xt_number(9)\n\
	 None<KeyPress>Return:	xt_goto_frame() xt_reset_number()";

static XtActionsRec window_actions[] = {
	{"PenRepaint",	      PenRepaint},
	{"xt_quit",           actionQuit},	
	{"xt_next",           actionNext},	
	{"xt_stretchy",       actionStretch},	
	{"xt_prev",           actionPrev},
        {"xt_run",            actionRun},
        {"xt_stop",           actionStop},
        {"xt_restart",        actionRestart},
	{"xt_faster",         actionFaster},
	{"xt_slower",         actionSlower},
	{"xt_number",         actionNumber},
	{"xt_reset_number",   actionNumReset},
	{"xt_goto_frame",     actionGoto},
	{"xt_print_coord",    actionCoord},
	{"xt_run_mode",       actionRunMode},
};

extern int xt_next_num;
