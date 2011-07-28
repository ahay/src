/*$

=head1 NAME

xtpen - VPLOT filter for the X toolkit

=head1 SYNOPSIS

xtpen  plot_file1 [plot_file2 ...] [options]

xtpen < plot_file [options]

graphics_program | xtpen [options]

Seplib versions:

Xtpen < Plot1.h [Plot2.h ...] [options]

SEP_graphics_program | Xtpen [options]

=head1 DESCRIPTION

B<pen> filters accept input in the B<vplot>
graphical metalanguage and cause a plot to appear on a device.
B<Xtpen> is the B<vplot> filter for the X toolkit.

If B<xtpen> is run without input or arguments it will
self-document. This is standard behavior for all
B<vplot> programs.

B<Xtpen> is an X toolkit program and
takes all the standard X-toolkit options,

E.g. -geometry 500x400 -font fixed

=head1 OPTIONS

Only xtpen-specific options are covered here.  For a list of all generic
B<pen>
options (which also apply to xtpen) do ``man vplot''. (If the documentation
you find is about Versatec plotters, that is an unrelated man page with
the same name.)

All B<xtpen>
options follow the standard SEPlib style, that is
on the command line in the format ``option=value'' with no
spaces on either side of the equal sign. You can also search for options in
a parameter file, by putting ``par=parfile'' on the command line (where
``parfile'' is the name of the parameter file).

=over 2

=item B<numcol=numcol>

   must be at least 8 and at most 256. It must be a power of 2.

   If left unspecified, B<xtpen> tries to request 128 free colors.  
   If it can't get that, then it tries for 64,
   then 32, and finally 16.
   (If it can't get even 16 it gives up and allocates a new color 
   table.) If B<xtpen>
   is only able to get a few color table entries (say, because 
   some X utility you run forgets to release color table entries 
   when it exits) then plots that need lots of colors may look 
   bad.

   To demand a specific number of colors,
   specify for example B<numcol=128> or B<numcol=256>
   on the command line, or put a line like
   B<``XTpen.numCol: 128''> in your X defaults file.
   However, note that if
   B<xtpen>
   is forced to allocate a new color table to get the requested 
   number of colors, the color table will only be correct 
   while the mouse is within the B<xtpen> window.

=item  B<buttons=yes> 

    If B<buttons=no,> the panel of 8 buttons at the top of the 
    plot will be ommitted.

=item  B<labels=yes> 

    If B<labels=no>, the display frame number and the inter-frame 
    pause time at the top of plot
    will be ommitted.

=item B<want_text=yes> 

   If B<want_text=no>,
   the message window will be ommitted.

=item B<message='text string'>

   This option, if present, displays the given text in 
   the message window.

=item B<see_progress=no> 

   If B<see_progress=y>,
   you can watch as each frame is plotted. This is slower, and 
   can cause problems if you attempt to move, obscure, or 
   iconify the plot while it is being drawn. You also may not 
   use a plotting window bigger than the actual screen if
   B<see_progress=y>.

   The X toolkit option B<-sync> can be used to dramatically 
   slow the speed of plotting.

=item B<stretchy=no> If B<stretchy=y>,

   xtpen will start in stretchy mode, with the default 10.24 
   high by 13.65 wide vplot plotting area anisotropically 
   stretched to fill the B<xtpen> window.  If B<stretchy=n>, 
   the default, B<xtpen> will follow the usual vplot-filter 
   practice; it will fit the largest 4 wide by 3 tall 
   rectangular plotting area it can within the B<xtpen>
   window, and will isotropically scale the vplot plotting 
   area into that.  (This assumes relative scaling; see the
   B<size> option below.)

=item B<images=yes pixmaps=no>

   These options are designed to allow fast redisplay of 
   multi-frame plots so that vplot files can be animated. 
   If they are both turned off the frames will be replotted 
   from scratch each time they are displayed.  Either 
   B<images=yes> or B<pixmaps=yes> will make B<xtpen> act 
   like a movie program, but in different ways.

   If B<images=yes>, B<xtpen>
   will copy the image created by plotting each frame and 
   save it in the client program B<(xtpen>).
   This will increase memory usage on the machine running
   B<xtpen>.
   If you have a fast connection to your X-server it will 
   make redisplay of frames faster. If you have a slow 
   connection e.g. telephone or appletalk it can actually 
   make replotting slower.

   If B<pixmaps=yes>, B<xtpen>
   will copy the image created by plotting each frame and 
   save it in the X-server.  Using this option,
   redisplay of frames will be very fast and the network 
   traffic will be very low; if your X-server is a workstation 
   with plenty of memory and swap space but your network 
   connection is slow, this is a useful option.

   However, if your X-server is an X-terminal or has 
   limited memory, this option
   may have undesirable effects on the response of your terminal.
   On many workstations memory allocated to the X server will 
   stay allocated and unavailable for general use until X windows, 
   not just B<xtpen>, is terminated.

=item B<interact='filename'>

   If this option is specified, at the end of plotting each frame
   B<xtpen>
   will give you the opporunity to click the cursor on points of 
   interest.  The vplot coordinates of each point will be written
   into the file 'filename'.

=item B<boxy=n> If B<boxy=y>, and B<interact>

   is set to a file name,
   a dialog box will pop up every time you click on the plot,
   giving you the opportunity to specify a label to go along with
   the picked coordinate. The output will be in a form
   suitable for input into the
   B<Box> program (used to annotate vplot plots).

=item B<ppi=(float) aspect=(float)>

   These options allow overriding of the true pixels per inch and 
   aspect ratio of the X display device. This can be useful if you 
   need to carefully align a plot with the pixels of a device, 
   regardless of the absolute size and shape of the device's 
   pixels.

=item B<size=relative>

   Normally when you plot something, you want the plot
   to scale to fill the page the same way it fills your 
   terminal screen.  This is B<size=relative>,
   the default behavior. If you wish vplot inches
   to really correspond to physical inches, you must 
   specify the generic B<pen> option B<size=absolute> instead.
   (Note the scale may be incorrect if the monitor is a 
   different size than what B<xtpen> thinks it is.)

   Note it is possible for the plot file to set the default 
   itself; in particular the vplot-to-vplot pen
   filter B<vppen>
   by default forces ``absolute style''
   for all of its output files.

=back

 

=head1 GENERIC OPTIONS

There are many many more generic B<pen>
options that B<xtpen>
also understands.
For example, to make everything fatter or thinner,
to rotate the plot, to rescale the plot,
to rescale only text, to use different fonts,
etc, etc, etc.

See the B<pen> manual page for a listing of all generic options.

Note B<xtpen>
defaults the universal option B<break>
to ``i'' (ignore) instead of the usual ``b'' (break)!

=head1  DISPLAY MODES

B<xtpen> has two display modes:

RUNNING MODE: Runs through all the frames in a loop;
active buttons are:

   QUIT : quits the program
   STOP : enter frame mode

FRAME MODE (pause=-1): Pauses after each frame; active buttons are:

   NEXT : go to next frame
   PREV : go to previous frame
   QUIT : exit the program
   RESTART : go to the first frame
   RUN  : enter running mode
   STRETCHY/RIGID : make plot fill the frame,
or preserve aspect ratio, respectively
   FORWARDS/BACKWARDS/BOTH-WAYS : change direction of frame flipping

Note that a backwards run will only show those frames already plotted.
It is advisable to run once through all the frames forwards.

=head1 X RESOURCES

The following actions are available for binding to keystrokes:

 xt_quit(): quit program
 xt_next(): next frame
 xt_prev(): prev frame
 xt_run(): run mode
 xt_stop(): frame mode
 xt_restart(): first frame
 xt_faster(): reduce pause between frames in run mode
 xt_slower(): increase pause between frames in run mode
 xt_stretchy(): toggle between stretchy and rigid modes
 xt_number(digit): enter a digit in the current number
 xt_reset_number(): reset the current number
 xt_goto_frame(): goto the frame number given by the current number
 xt_print_coord(): Print mouse coords in the file given by interact=

The default key bindings are:

 <Btn1Down>:   xt_print_coord()  \\n\\
 <KeyPress>n:  xt_stop() xt_reset_number() xt_next()  \\n\\
 <KeyPress>m:  xt_stop() xt_reset_number() xt_prev()  \\n\\
 <KeyPress>r:  xt_run()  \\n\\
 <KeyPress>q:  xt_quit()  \\n\\
 <KeyPress>.:  xt_stop()  \\n\\
 <KeyPress>f:  xt_faster()  \\n\\
 <KeyPress>s:  xt_slower()  \\n\\
 <KeyPress>t:  xt_stretchy()  \\n\\
 <KeyPress>0:  xt_number(0)  \\n\\
  ......                .......
 <KeyPress>9:  xt_number(9)  \\n\\
 <KeyPress>Return: xt_goto_frame() xt_reset_number() \\n\\
 <KeyPress>Escape: xt_reset_number()
Here is an example of overriding these in your ~/.Xdefaults file.
This binds the keypad number 1 to jump to the first frame:

 xtpen*translations: #override\\n\\
 <KeyPress>Q:    xt_quit() \\n\\
 <KeyPress>KP_1: xt_number(1) xt_goto_frame() xt_reset_number()

Here is an example of how to change a label:
 XTpen*Quit.label: QUIT

N.B.: using prev when you are at the first frame takes you to the last
frame plotted so far; this is not necessarily the last frame!
You can only jump to a frame if it has already been plotted once.
If you give an invalid frame number it will jump to the next frame.

Many parameters may have their defaults set in the X resource database.
Here are their equivalent names:

  option name  X-Resource name     Type
  ===========  ===============     ====
  stretchy      XTpen.stretchy     Boolean
  images        XTpen.useImages    Boolean
  pixmaps       XTpen.usePixmaps   Boolean
  buttons       XTpen.showButtons  Boolean
  labels        XTpen.showLabels   Boolean
  want_text     XTpen.showText     Boolean
  numcol        XTpen.numCol       int
  pause         XTpen.pause        int

E.g. If you want xtpen to come up in stretchy mode by default,
put this line in your Xdefaults file:

 XTpen.stretchy: True

=head1 USING THE ENTIRE SCREEN

The default XTpen app-defaults file that comes with the
B<xtpen>
source defines several options useful for photographing slides off the Sun
workstation screen, or for making demos. If these do not work for you,
you may need to define these in your personal app-defaults directory in
a file ``XTpen'':

 fullscreen.geometry: 1146x894-0-0
 slide.geometry: 1158x906-0-0
 XTpen.geometry: 800x700+300+125

If these are in your app defaults, then
B<xtpen -name fullscreen>
will plot using the largest window that will fit in the screen, and
B<xtpen -name slide>
will be even slightly larger, so that it is possible to push the window frame
off the screen leaving only the plot showing. (Note though that
if the plot is bigger than the screen, then
B<see_progress=y>
cannot be used.)

=head1  BUGS

If B<see_progress=y,>
each frame of the plot may attempt to redraw several times if the
plot is very small, or the first part of the plot may go missing if
the plot is very complex. This is probably an X toolkit bug, as it
appears to happen only with certain versions of X windows.

=head1 SEE ALSO

L<vppen>, L<pen>, L<pspen>, L<vplot>, L<vplottext>

=head1 COPYRIGHT

The Vplot source is copyrighted. Please read the copyright which can be
found in the accompanying
B<vplot> manual page.

=head1 AUTHOR

The device-independent code common to all vplot filters was primarily
written by Joe Dellinger, although many other people have made significant
contributions.
The device-dependent code of xtpen was almost entirely written by Dave Nichols.

=head1 CATEGORY

B<graphics/vplot/filters>

=cut

*/
const char *documentation[] = {
"NAME",
#ifdef SEP
"	Xtpen - SEPLIB vplot filter for X windows using the X Toolkit (Xt)",
#else
"	xtpen - vplot filter for X windows using the X Toolkit (Xt)",
#endif
"",
"SYNOPSIS",
#ifdef SEP
"	Xtpen [options] in=vplot-input file OR header file OR stdin",
#else
"	xtpen [options] [inputfiles]",
#endif
"",
"OPTIONS",
" This pen takes all the standard X-toolkit options",
" E.g. -geometry 500x400 -font fixed",
"",
" The pen has two display modes ",
"",
" RUNNING MODE: Runs through all the frames in a loop",
" Active buttons are: ",
"	QUIT : quits the program",
"	STOP : enter frame mode",
"",
" FRAME MODE (pause=-1): Pauses after each frame ",
" Active buttons are: ",
" 	NEXT : next frame",
"	PREV : previous frame",
"	QUIT : quits the program",
"	RESTART : go to the first frame",
"	RUN  : enter running mode",
"	STRETCHY/RIGID : make plot fill the frame or preserve aspect ratio",
"	FORWARDS/BACKWARDS/BOTH-WAYS : change direction of frame flipping",
"	Note that a backwards run will only show those frames already plotted",
"	It is advisable to run once through all the frames forwards.",
"",
" The following actions are available for binding to keystrokes;",
" xt_quit(): quit program   xt_next(): next frame   xt_prev(): prev frame",
" xt_run(): run mode        xt_stop(): frame mode   xt_restart(): first frame",
" xt_faster(): reduce pause between frames in run mode ",
" xt_slower(): increase pause between frames in run mode ",
" xt_stretchy(): toggle between stretchy and rigid modes",
" xt_number(digit): enter a digit in the current number ",
" xt_reset_number(): reset the current number ",
" xt_goto_frame(): goto the frame number given by the current number ", 
" xt_print_coord(): Print mouse coords in the file given by interact=",
"",
" The default key bindings are:",
"       <Btn1Down>:            xt_print_coord()  \\n\\ ",
"       <KeyPress>n:           xt_stop() xt_reset_number() xt_next()  \\n\\ ",
"       <KeyPress>m:           xt_stop() xt_reset_number() xt_prev()  \\n\\ ",
"       <KeyPress>r:           xt_run()  \\n\\ ",
"       <KeyPress>q:           xt_quit()  \\n\\ ",
"       <KeyPress>.:           xt_stop()  \\n\\ ",
"       <KeyPress>f:           xt_faster()  \\n\\ ",
"       <KeyPress>s:           xt_slower()  \\n\\ ",
"       <KeyPress>t:           xt_stretchy()  \\n\\ ",
"       <KeyPress>0:           xt_number(0)  \\n\\ ",
"        ......  		.......  ",
"       <KeyPress>9:           xt_number(9)  \\n\\ ",
"       <KeyPress>Return:      xt_goto_frame() xt_reset_number()  \\n\\ ",
"       <KeyPress>Escape:      xt_reset_number() ",
"",
"Here is an example of overriding these in your ~/.Xdefaults file",
" this binds the keypad number 1 to skip to the first frame ",
" xtpen*translations: #override\\n\\ ",
" <KeyPress>Q:       xt_quit() \\n\\ ",
" <KeyPress>KP_1:       xt_number(1) xt_goto_frame() xt_reset_number() ",
"",
" N.B using prev when you are at the first frame takes you to the last",
" frame plotted so far; this is not necessarily the last frame!",
" You can only jump to a frame if it has already been plotted once.",
" If you give an invalid frame number it will jump to the next frame.",
"",
" numcol=128 (must be power of 2)",
" buttons=yes (display a panel of buttons at the top of the plot)",
" labels=yes  (display frame number and inter-frame delay at the top of plot)",
" message=\" some text \" (displays the text in the message window)",
" want_text=yes (display a message window)",
" see_progress=no (force pen to show progress of each frame, slow)",
" stretchy=no (start pen in stretchy mode, vplot frame fills the window)",
"",
" The next two options are designed to allow fast redisplay of multi frame",
" plots so that vplot files can be animated. If they are both turned off the",
" frames will be replotted from scratch each time they are displayed.",
" using images=yes or pixmaps=yes will make the pen act like a movie",
" program.",
"",
" images=yes  Copy the image created by plotting each frame and save it in",
"     the client program (xtpen). This will increase memory usage in",
"     the machine that runs the pen command!",
"     If you have a fast connection to your X-server it will make redisplay",
"     of frames faster. If you have a slow connection e.g. phone or appletalk",
"     it may make replotting slower.",
" pixmaps=no  Copy the image created by plotting each frame and save it",
"     in the X-server. This will increase memory usage of the machine",
"     that displays the window!",
"     Redisplay of frames will be very fast and the network traffic is very",
"     low so this is a suitable option for slow connections.",
"     If your X-server is a workstation with plenty of memory and swap space",
"     then this option should be very useful.",
"     If your X-server is an X-terminal, or has limited memory, this option",
"     may have undesirable effects on the response of your terminal.",
"     On the other hand this the by far the fastest re-plotting mode!",
"",
" boxy=y interact=filename",
" Output the coordinates and labels suitable for input into the Box program",
" This option will popup a dialog every time you click on the plot",
"",
" interact=filename",
" This writes the vplot coordinates of every click in the file \"filename\".",
"",
" ppi=(float) aspect=(float) x_screen_info=no",
" These allow overriding of the true pixels per inch and aspect ratio",
" of the X display device. This can be useful if you need to carefully",
" align a plot with the pixels of a device, regardless of the absolute",
" size and shape of the device's pixels. Setting x_screen_info=yes outputs",
" the calculated default values.",
"",
"Many parameters may have their defaults set in the Xresource database",
"Here are the equivalent names:",
"  option name          X-Resource name         Type",
"  ===========          ===============         ====",
"  stretchy		XTpen.stretchy         Boolean",
"  images		XTpen.useImages        Boolean",
"  pixmaps		XTpen.usePixmaps       Boolean",
"  buttons		XTpen.showButtons      Boolean",
"  labels		XTpen.showLabels       Boolean",
"  want_text		XTpen.showText         Boolean",
"  numcol		XTpen.numCol           int",
"  pause		XTpen.pause            int",
"",
"E.g. If you want xtpen to come up in stretchy mode as a default",
"     put this line in your Xdefaults file:",
"XTpen.stretchy: True",
"",
#include "../include/gendoc.h"
"",
"Note xtpen defaults the universal option break to i (ignore), not b!",
"",
"SEE ALSO",
"	man pen"
};
int doclength = {sizeof documentation/sizeof documentation[0]};
