/*

=head1 NAME

pspen - VPLOT filter for Postscript

=head1 SYNOPSIS

pspen  plot_file1 [plot_file2 ...] [options] [> out.ps]

pspen < plot_file [options] [> out.ps]

graphics_program | pspen [options] [> out.ps]

Seplib versions:

Pspen < Plot1.h [Plot2.h ...] [options] [out=out.ps]

SEP_graphics_program | Pspen [options] [out=out.ps]

=head1 DESCRIPTION

All B<pen> filters accept input in the B<vplot>
graphical metalanguage, and either cause a
plot to appear on a device or write a file which if sent to that device will
create a plot.
B<Pspen> is the B<vplot>
filter for the Postscript page description language.
If called without
redirected output,
B<pspen>
creates a temporary file with a unique name,
spools the plot,
and removes the temporary file upon exiting.
If the output is redirected,
the Postscript output is written to standard out (out= for the Seplib
version) and no spooling is done.

If B<pspen>
is run without input, output, or arguments, it will
self-document. This is standard behavior for all
B<vplot> programs.


=head1  OPTIONS

Only pspen-specific options are covered here.  For a list of all generic
B<pen>
options, which also apply to pspen, do ``man vplot''. (If the documentation
you find is about Versatec plotters, that is an unrelated man page with
the same name.)

All
B<pspen>
options follow the standard SEPlib style, that is
on the command line in the format ``option=value'' with no
spaces on either side of the equal sign. You can also search for options in
a parameter file, by putting ``par=parfile'' on the command line (where
``parfile'' is the name of the parameter file).

The options below are listed in approximate order of relevance.

=over 2

=item wstype=

   If no ``workstation type'' is specified, a standard SPARC 
   laser printer is assumed. Other options currently known are
   B<wstype=oyo>, for an oyo electrostatic plotter, and
   B<wstype=tcpr>, for a Tektronix color plotter.
   (You can also achieve the same effect by running the 
   executable under the name ``oyopen'' or ``tcprpen'', 
   respectively.) Setting the workstation type determines the 
   default pixels per inch, plottable area, grey-scale 
   correction factors, etc.

=item printer=postscript

   If the output of pspen is not redirected into a file,
   by default B<pspen>
   spools the plot to the printer named ``postscript'' 
   on your system.  If that is not where you want the plot 
   to go, specify another printer name.

=item paper=

   Set the paper size. B<paper=letter>
   is the usual default (it can be set to something else
   at compile time, for example
   B<paper=a4> in Europe).
   Other possibilities are
   B<paper=legal,> B<paper=a3,> B<paper=a4,> and B<paper=a5.>

   It is also possible to set a custom paper size; for example
   B<paper=8.5x11i> is (almost) the same as B<paper=letter>.
   Usually the first number will be the shorter axis and the 
   second the longer; if you make either one 0 then
   B<pspen>
   will determine the other by assuming paper with a
   standard 4 by 3 aspect ratio.
   If the ``i'' at the end is left off inches are assumed. If you 
   prefer centimeters, use ``c''.

   By default
   B<pspen>
   assumes that the outermost .25 inches of the paper cannot 
   be plotted on.  If that is incorrect, you can override by 
   listing the two border widths
   after a comma after the units character.
   For example if there were no borders at all on the
   shorter axis of letter-size paper on your plotter you could 
   specify B<paper=8.5x11i,0.,.25>

   As a last resort for troublesome plotters it is possible to 
   specify the true minimum, true maximum, and theoretical 
   paper size for each axis of the paper explicitly; for a 
   SPARC printer B<paper=letter> is equivalent to 
   B<paper=.27,8.23,8.5x.19,10.88,11.>
   (These numbers were determined empirically.)

=item size=relative

   Normally when you plot something, you want the plot
   to scale to fill the page the same way it fills your 
   terminal screen.  This is B<size=relative,>
   the default behavior. If you wish vplot inches
   to really correspond to physical inches, regardless of 
   the size of the paper, you must specify the generic
   B<pen> option B<size=absolute> instead.

   Since vplot's standard ``plotting window'' is
   10.24 by 13.65 inches, for 8.5 by 11 paper these two 
   scaling options will not produce very different results.
   For much larger or smaller paper than vplot's plotting 
   window, however, the results will be very different.
   It is also possible for the plot file to set the default 
   itself; in particular the vplot-to-vplot pen filter 
   B<vppen> by default forces ``absolute style''
   for all of its output files.

=item label='string'

   By default, plots are labeled in the lower right corner
   with the username of the person making the plot and the time 
   and date.  This can be replaced with something else by using 
   the B<label> option, for example B<label='Dellinger/Muir 
   Figure 3a'.> (If the string includes spaces, be sure to 
   enclose the string in quotes.)

   The label may obscure part of the plot, as it whites out
   a rectangle in the area it will occupy before writing text 
   in it.  B<label=''> or B<label=no>
   turns off the label and the white-out rectangle.


=item color=n force=n invras=y 

   If B<color=n,>
   the output postscript file assumes a monochrome printer.
   Any non-background color will be plotted as black.
   Raster will be converted to grey-scale and then dithered to 
   simulate a grey-scale.
   Polygons will either be filled with solid black, solid white, 
   or perhaps some dither pattern (depending on how the polygon 
   is defined in the vplot input file, and the particular 
   implementation of postscript on the printer).

   If B<color=n,> the generic pen option B<invras=n>
   will invert just the raster grey-scale.  See the B<pen>
   man page for more details.

   If B<color=y,> B<pspen>
   produces a color postscript output.
   Noncolor printers will produce
   grey-scale output; B<grey=y> is a synonym for B<color=y>
   for backwards compatibility.

   Since on screen devices the background is black
   while on printers the background is white, by default
   B<pspen>
   complements all non-raster  colors to preserve their relative 
   intensity (white becomes black, blue becomes yellow, red 
   becomes cyan, etc).  If this is not desired,
   B<force=y> forces the non-raster color to come out true
   (although the background will still be white).
   This option has no effect unless
   B<color=y.>


   Note that if B<color=y.>
   colors inside raster blocks are always shown true, regardless
   of the setting of B<force!>

   Generally if you are not resetting the background color then 
   the default, B<force=n,> will produce a reasonable result. 
   If you are resetting the background (say, to dark blue) then 
   you might consider using B<force=y> and explicitly beginning 
   your plot by drawing a background rectangle of the desired 
   color.


=item strip=no

   If B<strip=y,>
   then the output Postscript file has no label or document 
   control language prepended/appended.
   This makes the output postscript file suitable for inclusion 
   into other postscript documents.
   Using the shell ``pstexpen'' you can also generate a 
   Postscript file with the bounding box information at the 
   start that many editors, word processors, etc, expect to 
   find when they import postscript.  B<tex=y> is a synonym 
   for B<strip=y> provided for backwards compatibility.

=item hold=no

   If B<hold=y,> the printer will expect manually fed paper 
   for this plot.

=item ncopies=1

   This parameter is used to save time and effort by printing 
   multiple copies of the plot.

=item ppi=300. 

   B<Pspen>
   chooses a suitable ``pixels per inch'' depending on the 
   workstation type.  (Currently 300 dots per inch for all 3 
   specifically supported devices.) If this is incorrect, 
   override the default pixels per inch from the command line.

=item greyc=-0.5 pixc=0.6

   Grey scale reproduction on a hardcopy device is quite 
   different from that on a display device. The transition 
   from black to white occurs more abruptly on a display 
   device, leaving both ends of the grey scale clipped at 
   black or white. This nonlinearity in the perceived grey 
   scale can be simulated on the printer output using the
   B<greyc> parameter.
   Even after taking the contrast correction
   into account, the neutral-grey parts of a grey scale image 
   usually come out too dark when plotted on hardcopy compared 
   to a good graphics screen. This darkening is
   caused by the overlapping and bleeding of nearby black dots, 
   distorting and darkening the printer dithering pattern.
   The B<pixc> parameter compensates for this pixel overlap by 
   shifting the darker half of the color scale
   towards white.

   Theoretically printer engines are supposed to take
   all these factors into account, but we have found that
   empirically B<greyc=-0.5 pixc=0.6>
   seems to give truer grey scale reproduction on most 
   postscript printers.  B<pixc=1. greyc=1.> disables this 
   correction.  Other printer engines, notably ``oyo'', may 
   require very different values (see the ``wstype'' option).
   Further experimentation may
   be required to determine suitable values for your particular 
   printer engines.

   In particular, when preparing SEG extended abstracts 
   fiddle with pixc and greyc until the raster looks very 
   washed-out and light; otherwise, when shrunk and printed 
   in the abstract volume your grey scale plots will
   become much too dark.

   See the B<vplotraster> manual page
   for more details.

=item harddither=y dither=4

   Dithering is a means of representing a continuous-tone 
   grey image on a hardcopy device. Usually it is best to 
   use the postscript printer's built-in hardware dithering
   B<(harddither=y),>
   but if that is poorly implemented you may want to use vplot's
   software dithering instead
   B<(harddither=n).>

   If vplot does the dithering there are several methods 
   to choose from.  The default dithering method is the
   minimized average error algorithm
   B<(dither=3).>
   This dithering method provides the highest resolution
   but does not photocopy well.
   You might prefer the digital halftoning algorithm
   B<(dither=4)> instead.  (This is similar to the built-in 
   algorithm used by most postscript 
   printers.) See the B<Vplotraster> manual page for a 
   complete discussion of dithering options.

=item dumbfat=n 

   If plot will emulate drawing fat lines and dashed lines 
   in software.  (Some postscript implementations are buggy, 
   and barf on plots with fat lines in them.)

=back

=head1  OTHER OPTIONS

There are many many more generic
B<pen> options that B<pspen> also understands.  See the
B<pen> manual page for a listing of all generic options.
Various plot manipulations can be performed via
B<vppen>.  See the B<vppen> manual page for more details.

Some generic
B<pen> options you may find useful:

=over 2

=item fat=1 

   Make all vectors a little thicker.

=item fatmult=0. 

   Make all vectors thin.

=item scale=2.  

   Make the plot twice as big.

=item xscale=2. 

   Make the plot twice as wide.

=item yscale=2. 

   Make the plot twice as tall.

=item txsquare=y 

   Force text to not be stretched, despite B<xscale> and B<yscale.>

=item hshift=2

   Shift the plot to the right 2 inches.

=item vshift=2. 

   Shift the plot up 2 inches.

=item xcenter=0. ycenter=0.

   Place the vplot origin at the center of the plotting area.

=item txscale=2. 

   Make all text twice as big.

=back

=head1 ENVIRONMENT VARIABLES

=over 2

=item DEFAULT_PAPER_SIZE

   What size paper to use; equivalent to the ``paper'' 
   command line option.

=item WSTYPE

   Workstation type; equivalent to the ``wstype'' command 
   line option.

=item PSPRINTER

   Printer name; equivalent to the ``printer'' command line 
   option.

=item VPLOTFONTDIR

   The directory where special-purpose runtime-loaded vplot fonts 
   can be found.  Normally you won't need these fonts unless you 
   want to (for example) label your plots in Russian, or use the 
   occasional Greek letter, etc.

=item VPLOTSPOOLDIR

   Where to store temporary files. Normally ``/tmp''. If that 
   is not a good place, specify another.

=back

=head1 COMPILE-TIME DEFINES

=over 2

=item NOCOMMENTS 

   If you are used to having pspen spool plots silently,
   define NOCOMMENTS at compile time.

=item DEFAULT_PAPER_SIZE

   If you prefer a different default paper size, define this 
   to be the appropriate string (be sure to include double 
   quotes).

=item DEFAULT_PAPER_UNITS

   If you prefer default units of centimeters instead of 
   inches, define this to be `c' (with single quotes).

=back

=head1 BUGS

Pspen really should wait to pass through line dashing information until
it sees it's actually going to get used.

Some postscript interpreters don't interpret all postscript commands
(Ghostscript if you try to use a dithered polygon).
Some postscript interpreters have a nervous breakdown if they are
given a plot with too many vectors in it (Oyo and Versatec bed plotters).
Some postscript printers have memory leaks and slowly slow down
plot by plot until they are rebooted (some Apple laserwriters).
Some postscript interpreters don't do line fattening properly (some
SPARC interpreters; this is the cause of ``fanged S'' syndrome).

=head1 SEE ALSO

L<vppen>, L<pen>, L<vplot>, L<vplotraster>

=head1  COPYRIGHT

The Vplot source is copyrighted. Please read the copyright which can be
found in the accompanying
B<vplot> manual page.

=head1 AUTHOR
The device-independent code common to all vplot filters was primarily
written by Joe Dellinger, although many other people have made significant
contributions.
Pspen was originally written by Wes Monroe. Steve Cole has added several
new routines and extensively overhauled most of the old ones.
Over the years many other people have made minor changes.

=head1 CATEGORY

B<graphics/vplot/filters>

=cut

*/
const char *documentation[] = {
" ",
"NAME",
#ifdef SEP
"    Pspen - SEPlib vplot filter for Postscript",
#else
"    pspen - vplot filter for Postscript",
#endif
" ",
"    output is in PostScript language; if not redirected, it is sent to",
"    lpr -Ppostscript   (override with $PSPRINTER environment variable.)",
" ",
"SYNOPSIS",
#ifdef SEP
"    Pspen [options] in=vplot-inputfile OR headerfile on standard in",
#else 
"    pspen [options] [inputfiles]",
#endif
" ",
"OPTIONS",
#ifdef SEP
"    Options unique to Pspen:",
#else
"    Options unique to pspen:",
#endif
"   printer=postscript   (what printer to send it to)",
"   paper=letter   If paper=legal, 8.5x14 inch paper in use.",
"                     paper=a4 is also recognized.",
"                     For custom paper size do for example paper=8.5x11i",
"   label=string   Label pages with string.  If string=no, no label",
"                     default label is user name and date",
"   strip=no       If strip=y, then make an unrotated PostScript output file",
"                     without document control language",
"   hold=no        Tells the printer to not print the job until you",
"                     add paper through the manual feed slot.",
"   ncopies=n      Requests multiple copies (works only if strip=no)",
"   color=n        If y, then use COLOR postscript! (Monochrome printers will",
"                     use some fixed number of dither patterns for grey levels)",
"   force=n        If y, don't replace colors with their complements",
"",
" See pspen man page for obscure pspen options and environment variables.",
#include "../include/gendoc.h"
"FILES",
"    scratch file /tmp/PSPEN_XXXXXX, which it deletes",
"SEE ALSO",
"    man pspen pen vplot vplotraster"
};
int	doclength = { sizeof documentation/sizeof documentation[0] };
