/*
=head1 NAME

vppen - VPLOT filter for VPLOT

=head1 SYNOPSIS

vppen  plot_file1 [plot_file2 ...] [options] > vplot_out

Seplib version: Vppen < Plot1.h [Plot2.h ...] [options] > Vplot_out.h

=head1  DESCRIPTION

All B<pen> filters accept input in the B<vplot>
graphical metalanguage. B<Vppen>
is the filter for the ``virtual vplot device''.
Both its input and output are vplot. Why is this useful?
The full set of pen options allows you
to rotate, combine, clip, scale, shift, etc, plots. At some point you may
make something on your screen using these options which
you would really like to have as a single
file. How to do this? Simply use the same set of options with vppen.
Vppen will draw the same picture, but will save the picture as a vplot file.
You can then rotate, combine, clip, etc this new plot.

Vppen also has options to do useful things such as centering, making
arrays of plots, and finding plot statistics.


=head1 VIRTUAL DEVICE DESCRIPTION

The virtual vplot device has a screen 10.24 inches (STANDARD_HEIGHT in <vplot.h>)
high, with a height/width ratio of 3/4 (SCREEN_RATIO in <vplot.h>). It has
600 pixels to the inch (RPERIN), and has 256 colors, all settable.
It knows how to do its own clipping and has ``hardware text'', etc.
This is the device that vppen ``plots'' on.

Ideally vplot coming into
vppen with no options set should be passed on through completely unchanged.
This is true (except for a few exceptions like polygon fill patterns and such)
if both the input and output are absolute style. (If you don't know what
terms like ``absolute style'' mean, read the other vplot documentation.)
For this reason, vppen by default uses ``size=absolute'' and outputs absolute style
output. This is convenient if you are repeatedly sending files through
vppen, but you may want to override it otherwise.

=head1 GENERIC OPTIONS

For a list of all generic
pen options, which also apply to vppen, do ``man pen''.
Generic options that you may find particularly useful are:

=over 2

=item B<erase>

   The erase=once option is useful when combining several 
   plots into one.  (Otherwise you're likely to have the 
   entire screen cleared before each plot.)

=item B<xshift, yshift>

   These options are useful for positioning a plot.

=item B<scale>

   The scale options are useful in getting a plot the right size.

=item B<xcenter, ycenter, interact>

   When you're going crazy trying to get something centered 
   in an oddly shaped screen made by ``gridnum'', this 
   comes in very handy.

=item B<xwmin, xwmax, ywmin, ywmax>

   Useful for clipping out just the part of the plot you want.

=item B<rotate>

   Useful for rotating a plot, usually by 90 degrees.


Remember that ALL generic pen options apply, not just the ones listed
above.
Only vppen-specific options will be covered in the next section.

=back

 

=head1  VPPEN-SPECIFIC OPTIONS

=over 2

=item B<gridnum=0 gridsize= grid=-1>

   These commands are provided as a quick way to make a grid 
   of plots.  Gridnum=X,Y divides the output screen into 
   X rectangles horizontally and Y rectangles vertically. 
   If Y is not specified, then it is set equal to X. If X is 
   zero, then gridding is not done.  Gridsize=X,Y sets the X 
   and Y dimensions of each rectangle in inches.
   If the gridsize is not set, then the array of rectangles 
   will be sized to fill the screen.
   (The ``screen'' has an X of 13.65 inches and a Y of 10.24 
   inches.) Grid sets the fatness of a white border drawn 
   around each rectangle.  If it is negative, no border will 
   be drawn.  Each rectangle is its own miniature device 
   screen. Clipping is set to the edges of this screen just 
   as it is for the full-sized screen.  All generic commands 
   such as rotate, xshift, etc, will apply to each rectangle 
   and NOT to the plot as a whole.  Each rectangle contains 
   one (or zero) plot. When an erase command is received, 
   nothing will be erased but instead plotting will begin 
   in the next rectangle, and a ``break'' command will be 
   output.  The rectangles are filled left to right and then 
   top to bottom. It is not an error to have more plots that 
   you have room for, the rectangles will march right off the 
   bottom of the page, but they will still be there.

   If you have less plots than you have room for some of the 
   rectangles will be empty. Normally the generic pen option 
   ``size'' defaults to ``absolute'' in vppen, as it does on 
   most hardcopy devices. This doesn't work very well
   if your screen is a few inches on a side, so if the gridnum 
   option is used ``size'' will instead default to the normal 
   screen-device default ``relative''.  For similar reasons, 
   the gridnum option will also change the default options to 
   ``big=n'' and ``vpstyle=n''.  You may find it useful to 
   use size=relative when plotting a gridded vplot file.

=item B<big=y>

   If big=y, the output screen is expanded to the entire 
   range of possible vplot coordinates, -54.6 to +54.6 
   inches (VP_MAX in <vplot.h>) in each direction. A 
   coordinate shift is thrown in to move the origin back 
   to (0,0) where it belongs. This has some good effects 
   and some bad effects.  The chief good effect is that 
   nothing is clipped at the edge of the screen, since the 
   ``screen'' contains all possible vplot coordinates.  
   (It is still possible for things to be clipped, but only 
   at the edge of the vplot coordinate system.) The chief bad 
   effect is that rotated style or relative size objects
   will be positioned at the top of the expanded screen, which 
   is to say 54 inches away from the origin and right at the 
   edge of the largest coordinates vplot can handle.
   Thus big=y with rotated style or relative size
   input is an all around bad idea; the
   solution for rotated style plots is
   to make one pass through vppen with big=n to turn it
   into non-rotated style. There is no reason to specify 
   size=relative with big=y, unless you want a very big plot 
   indeed!

=item B<stat=n>

   If stat=y, no vplot will be written to standard out. 
   Instead vppen will find the largest and smallest coordinates 
   used in each input file and print out a summary of the height, 
   width, etc, information. Statistics are also kept for all the 
   input plots together as a whole.

=item B<align=uu>

   This option is used to left, right, top, bottom, or center 
   justify a plot.  The format is align=xy, where ``x'' controls 
   the left-right justification and ``y'' controls the up-down 
   justification. ``X'' is one of:

   l for left justified
   r for right justified
   c for centered
   u for unaligned, (no shifting done)

   and ``Y'' is one of:

   b for bottom justified
   t for top justified
   c for centered
   u for unaligned, (no shifting done)

The align point is set to have coordinate (0,0). Note that points shifted
into negative coordinates are still there, they just may be off the
screen. (Use the big=y option to avoid clipping problems.)
The ``xcenter=0, ycenter=0'' option is very handy to use when
plotting ``aligned'' files.

=item B< xsize= ysize=>

   These options allow you to scale your plot to fit within a
   rectangular area of the desired size. If you specify both xsize
   and ysize, the size of the area will be xsize inches wide by 
   ysize inches high. Such scaling can change the aspect ratio of 
   the plot.  If you specify EITHER xsize or ysize, vppen picks 
   the other so that the original aspect ratio is maintained.

=item B<vpstyle=y>

   Normally vppen inserts a ``set style absolute'' vplot 
   command at the beginning of every plot frame. If vpstyle=n, 
   it does not do this.

=item B<dumb=n>

   If dumb=y, then all output will be digested down to the 
   bare essentials - color changes, erases, moves, draws, 
   and nothing else.

=item B<blast=y bit=0>

   If blast=n, then raster output will be compacted as much 
   as possible as it is being written. This is slower, but 
   the resulting file may be substantially smaller. If 
   bit=integer > 0, then bit raster will be used, with the 
   integer giving the color number of the ``on'' pixels.

=item B<outN= outN+=>

   This option is used to redirect a given frame of the output 
   into a separate file. ``N'' is the frame number (the first 
   frame being number 0).  If for example ``out0=frame0'' were 
   found on the vppen command line, then vppen would attempt 
   to write the first frame's worth of vplot output into the 
   file ``frame0'' instead of the usual place. Optionally there
   can be a ``+'' after the frame number; in that case vppen 
   opens the output file for appending instead of attempting 
   to overwrite it.  Alternatively, use literally ``outN'' or
   ``outN+''; in that case the filename will be used in a
   printf-style format with the integer frame number as an
   argument.  For example, ``outN=file%d'' would put frame 0
   into ``file0'', frame 1 into ``file1'', etc. (An individually
   specified file name, with N an integer, will override.)
   Note that Since color table values are ``global'' from one
   plot frame to another, vppen explicitly writes out all the
   color table entries currently in use whenever the vplot output
   stream is redirected. Each frame will thus come out with the
   correct color table even when plotted individually.

=back

 

=head1 COMMENTS

Remember setting ``big=y'' and ``size=relative'' at the same time
is a very bad idea!

Some options (stat, align) will not work with piped input, as they involve
making repeated passes through the input vplot files.

Vppen should refuse to dump raw binary vplot to your screen.

Vppen outputs an initial erase only if an initial erase is forced from the
command line (erase=yes (the default) or erase=once).

Set the text font, precision, etc, that you want on your first pass
through vppen, because it will be hardwired after that.

Some vplot commands as yet have no corresponding command in libvplot, and
so vppen currently eats these. Provisions are made in the code for the
day when these commands will exist; they just have to be uncommented out.

Vppen has four uses. One is as a sort of cheap vplot editor. The second
is as a filter for digesting complicated things into the bare essential
moves and draws. The third is that it provides a library of routines
that can be linked into other vplot filters (perhaps an interactive vplot
editor), so that they can easily ``dump a hardcopy of the screen''.
(The code has been carefully written with this in mind.)
Lastly, it provides an example of how to do several tricky things:
First, how to support a device that can do
EVERYTHING in ``hardware'', and so needs to have the highest possible
level of support. Second, how to do multiple passes through the vplot
input (again something paving the way for a vplot editor).
Third, how to find the length of
a text string without actually drawing it.

=head1 SEE ALSO

L<pen>, L<libvplot>, L<vplot>

=head1 COPYRIGHT
The Vplot source is copyrighted. Please read the copyright which can be
found in the accompanying Vplot manual page.


=head1  AUTHOR
Joe Dellinger

=head1  BUGS
There are still rotated style plots in common use at SEP, which have
problems with the default big=y. There's not much I can do about this
one, either, and still be backwards compatible. See the discussion
under the ``big'' option to find out how to get around this problem.

The ``stat=y'' option
doesn't work very well with the SEPlib version of Vppen, as the statistics
go to the output data file, not to the screen. You can force the output
to the screen by specifying ``out=stdout'' or ``out=/dev/tty''.

It is not clear how options such as ``stat'', ``align'', and ``gridnum''
behave when used in combination.

Absolute sizing can be a pain. You may just want to alias ``vppen'' to
``vppen vpstyle=n''.

=head1 CATEGORY

B<graphics/vplot/filters>

=cut

*/

char *documentation[] = {
" ",
"NAME",
#ifdef SEP
"    Vppen - SEPlib vplot filter for the virtual vplot device",
#else
"    vppen - vplot filter for the virtual vplot device",
#endif
" ",
"Although it is perhaps not obvious, this program can be used to",
"\"Capture the screen\". Ie, you play with Pen options until you",
"get something you like, and then you can use those options with",
"this program to make a new vplot file that without any options",
"will draw the same thing.",
" ",
"OPTIONS",
"    Defaults: dumb=n stat=n big=y align=uu vpstyle=y blast=y bit=0",
"	       xsize=,ysize=,gridnum=0,0 gridsize=xwidth,ywidth grid=-1",
"	       out0= out0+= out1= out1+= out2= [etc]",
"dumb=y causes output to only be vectors, erases, and color changes.",
"stat=y causes plot statistics to be printed to stdout instead of vplot.",
"big=y expands the size of the device's screen (and hence outermost",
"clipping window) to nearly infinity (bad for rotated style!).",
"align=xy aligns plot:",
"x is one of l, r, c, u for left, right, center, unaligned",
"y is one of b, t, c, u for bottom, top, center, unaligned.",
"In all cases the given point is aligned to have coordinate zero.",
"xsize and ysize allow you to scale the vplot image to fit within",
"a given size rectangle. If you specify both the image may be",
"distorted by this scaling. If you specify only one of the two,",
"vppen selects the other so that the aspect ratio of the plot is",
"kept constant",
"vpstyle=n omits declaring absolute style in the output file.",
"(The generic pen option \"size\" defaults to absolute for input.)",
"blast is as in the libvplot raster documentation. ",
"if bit > 0, then bit raster is used with bit the color.",
"gridnum=numx,numy grids the screen, each part has gridsize=xwidth,ywidth",
"numy defaults to numx. [xy]size default to [xy]screensize / num[xy].",
"grid=N turns on drawing a grid, with fatness N.",
"",
"outN[+]=(filename) allows frame \"N\" of the output to be redirected",
"into a separate file. (The first frame of the plot is set by \"out0\".)",
"If there is a \"+\" after the frame number, the file will be opened for",
"appending instead of being overwritten. If you want all frames to be",
"redirected then use a literal \"outN\" and include a %d in the filename",
"template on the right where you want the frame number to go.",
"",
"Some combinations of options are not allowed.",
"You may not redirect or pipe the input if either the stat or align option",
"is used. ALL GENERIC PEN OPTIONS ARE APPLICABLE.",
};
int	doclength = { sizeof documentation/sizeof documentation[0] };
