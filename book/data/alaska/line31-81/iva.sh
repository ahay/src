#!/bin/sh
# File: iva.sh

##############################################################
# Credits:
# 2004 Hale, Cohen, with Stockwell modifications 2004.
#      In su distribution directory:
#      ~/su/src/demos/Velocity_Analysis/Traditional/Velan
# 2005 Seismic Processing with Seismic Un*x 
#      Forel, Benz, Pennington, 2005
#      script iva.sh  (section 7.6.7.3) and 
#      velanQC.sh (section 8.2.2.2)
# 2011 Schleicher, offerred to David Forel for new edition
#      of "Seismic Processing with Seismic Un*x"
##############################################################

# Set messages on
#set -x

#================================================
# USER AREA -- SUPPLY VALUES
#------------------------------------------------
# CMPs for analysis

eval $1
echo parm1 is $1
echo  cmp1 and numCMPs $cmp1 $numCMPs
#------------------------------------------------
# File names

indata=$2
outpicks=$3  # ASCII file
inpicks=$4

#######################################################3
#######################################################3
#------------------------------------------------
# display choices

myperc=98       # perc value for plot
SUXWIGB_OR_XIMAGE=suxwigb
SUXWIGB_OR_XIMAGE=suximage

# size of the display windows
HBOX=700    # originally 450
WBOXCVS=300 # originally 300
WBOXCMP=200 # originally 300
WBOXVELAN=500 # originally 300

XBOXVELAN=10 # kls 1350 puts it on second screen

# these parameters work nicely to put plots on
# my 2nd screen. I turn them on by changing 
# the next line to:    if [ 1 ] 
if [ 2 -eq 1 ]  
then
    echo "******** in if ***"
    HBOX=1000    # originally 450
    WBOXCVS=300 # originally 300
    WBOXCMP=200 # originally 300
    WBOXVELAN=800 # originally 300
    XBOXVELAN=1350
fi

XBOXCVS=`expr $XBOXVELAN + $WBOXVELAN`
XBOXCMP=`expr $XBOXCVS + $WBOXCVS`
XBOXNMOCMP=`expr $XBOXCMP + $WBOXCMP`
#------------------------------------------------
# Processing variables

# Semblance variables
nvs=101   # number of velocities
dvs=$6   # velocity interval
fvs=$5    # first velocity

# Compute last semblance (velan) velocity
lvs=`echo "$fvs + (( $nvs - 1 ) * $dvs )" | bc -l`

# CVS variables
fc=$fvs # first CVS velocity
lc=$lvs # last CVS velocity - now same at the last semblance velocity
nc=11   # number of CVS velocities (panels)
SPAN=11   # ODD number of CMPs to stack into central CVS

#================================================

# HOW SEMBLANCE (VELAN) VELOCITIES ARE COMPUTED

# Last Vel =  fvs + (( nvs-1 ) * dvs ) = lvs
#     5000 =  500 + ((  99-1 ) * 45  )
#     3900 = 1200 + (( 100-1 ) * 27  )

#------------------------------------------------

# HOW CVS VELOCITIES ARE COMPUTED

# dc = CVS velocity increment
# dc = ( last CVS vel - first CVS vel ) / ( # CVS - 1 )
# m = CVS plot trace spacing (m = d2, vel units)
# m = ( last CVS vel - first CVS vel ) / ( ( # CVS - 1 ) * SPAN )

# j=1
# while [ j le nc ]
# do
#   vel =  fc  + { [(   lc - fc   ) / ( nc-1 )] * ( j-1) }
#   j = j + 1
# done
# EXAMPLE:
#   vel = 1200 + ( (( 3900 - 1200 ) / ( 10-1 )) * ( 1-1) )
#   vel = 1200 + ( (( 3900 - 1200 ) / ( 10-1 )) * ( 2-1) )
#   ...
#   vel = 1200 + ( (( 3900 - 1200 ) / ( 10-1 )) * (11-1) )
#================================================

# FILE DESCRIPTIONS

# spanpanel.$picknow.su = binary temp file for input CVS gathers
# cvs.$picknow.su = binary temp file for output CVS traces
# nmopanel.$picknow.su = binary temp file for NMO (flattened) section
# panel.$picknow.su = current CMP windowed from line of CMPs
# spanpanel.$picknow.su = span of CMPs windowed to make cvs
# picks.$picknow = current CMP picks arranged as "t1 v1"
#                                                "t2 v2"
#                                                 etc.
# par.# (# is a sequential index number; 1, 2, etc.)
#      = current CMP picks arranged as
#        "tnmo=t1,t2,t3,...
#        "vnmo=v1,v2,v3,...
# par.0 = file "par.cmp" re-arranged as
#         "cdp=#,#,#,etc."  NOTE: # in this line is picked CMP
#         "#=1,2,3,etc."    NOTE: # in this line is "#"
# outpicks = concatenation of par.0 and all par.# files.

#================================================


echo "  *** INTERACTIVE VELOCITY ANALYSIS ***"


#------------------------------------------------
#kls Remove old files.  Open new files
#rm -f panel.*.su picks.* par.* tmp*

echo "save the old outpicks file with date in the name:"
suffix=`date|sed "s/ /_/g"`
echo  "cp -p $outpicks $outpicks.$suffix"
cp -p $outpicks $outpicks.$suffix

#make a set of picks.* files that will be used to plot on the
#velan and the cvs plots, and apply nmo to gather plot

if [ -s $inpicks ]
then
 tvnmoqc mode=2 prefix=picks par=$inpicks
fi

#------------------------------------------------
# Get ns, dt, first time from seismic file
nt=`sugethw ns < $indata | sed 1q | sed 's/.*ns=//'`
dt=`sugethw dt < $indata | sed 1q | sed 's/.*dt=//'`
delrt=`sugethw delrt < $indata | sed 1q | sed 's/.*delrt=//'`

# Convert dt from header value in microseconds
# to seconds for velocity profile plot
dt=`echo "scale=6; $dt / 1000000 " | bc -l`

# If "delrt", use it; else use zero
tstart=`echo "scale=6; ${delrt} / 1000" | bc -l`

#------------------------------------------------
# BEGIN IVA LOOP
#------------------------------------------------

i=1

while [ $i -le $numCMPs ]
do
 # set variable $picknow to current CMP
 eval pickprev=\$cmp`expr i - 1`
 eval picknow=\$cmp$i
 eval picknext=\$cmp`expr i + 1`

 # make a file so the while loop will run the first time
 echo "just some junk" > newpicks.$picknow
 while [ -s "newpicks.$picknow" ] 
 do
  # work this location until the user makes no picks
  # on the velan
  
  if [ -s picks.$picknow ] ; then
    echo "Location CMP $picknow has no picks."
  fi

  #------------------------------------------------
  # Plot CMP (right)
  #------------------------------------------------

  suwind < $indata \
           key=cdp min=$picknow max=$picknow > panel.$picknow.su

  $SUXWIGB_OR_XIMAGE < panel.$picknow.su \
      xbox=$XBOXCMP ybox=0 wbox=$WBOXCMP hbox=$HBOX \
      title="CMP gather $picknow" \
      label1=" Time (s)" label2="Offset (m)" key=offset \
      perc=$myperc verbose=0 &

  #------------------------------------------------
  # Constant Velocity Stacks (CVS) (middle-left)
  # Make CVS plot for first pick effort.
  # If re-picking t-v values, do not make this plot.
  #------------------------------------------------

  # truncate SPAN to odd number less then or equal to SPAN
  HALF_SPAN=`expr $SPAN / 2`
  SPAN=`expr $HALF_SPAN "*" 2 + 1`

  # Select CMPs $picknow +/- HALF_SPAN. 
  # Write to  spanpanel.$picknow.su
  CMPMIN=`expr $picknow - $HALF_SPAN`
  CMPMAX=`expr $picknow + $HALF_SPAN`
  suwind < $indata key=cdp min=$CMPMIN max=$CMPMAX \
  > spanpanel.$picknow.su

  # Calculate CVS velocity increment
  # dc = ( last CVS vel - first CVS vel ) / ( # CVS - 1 )
  dc=`echo "( $lc - $fc ) / ( $nc - 1 )" | bc -l`

  # Calculate trace spacing for CVS plot (m = d2, vel units)
  # m = ( last CVS vel - first CVS vel ) / ( ( # CVS - 1 ) * SPAN )
  m=`echo "( $lc - $fc ) / ( ( $nc - 1 ) * $SPAN )" | bc -l`
  if [ ! -s cvs.$picknow.su  ] ; then
    # CVS ve locity loop
    rm cvs.$picknow.su
    j=1
    while [ $j -le $nc ]
    do
      vel=`echo "$fc + $dc * ( $j - 1 )" | bc -l`

      # uncomment to print CVS velocities to screen
      ##    echo " vel = $vel"

      sunmo < spanpanel.$picknow.su vnmo=$vel |
      sustack >> cvs.$picknow.su

      j=`expr $j + 1`
    done
  fi

  # Compute lowest velocity for annotating CVS plot
  # loV = first CVS velocity - HALF_SPAN * vel inc
  loV=`echo "$fc - $HALF_SPAN * $m" | bc -l`

  #------------------------------------------------
  # Picking instructions
  #------------------------------------------------
  echo " "
  echo "Preparing CMP $i of $numCMPs for Picking "
  echo "Location is CMP $picknow. CVS CMPs = $CMPMIN,$CMPMAX"
  echo " "
  echo "  Use the semblance plot to pick (t,v) pairs."
  echo "  Type \"s\" when the mouse pointer is where you want a pick."
  echo "  Be sure your picks increase in time."
  echo "  To control velocity interpolation, pick a first value"
  echo "    near zero time and a last value near the last time."
  echo "  Type \"q\" in the semblance plot when you finish picking."
  echo " "
  echo " If there are no picks (using \"s\") before you quit (using"
  echo " \"q\" in the semblance plot, picking will continue at the"
  echo " next CMP."  
  #------------------------------------------------
  # Plot semblance (velan) (left)
  #------------------------------------------------

  # if there is a non-zero length picks.$picknow file, plotit
  # kls add logic for pickprev, picknow, picknext
    if [ -s picks.$picknow ]
    then

      #---  ---  ---  ---  ---  ---  ---  ---  ---  ---
      # Get the number of picks (number of lines) in picks.$picknow
      #   Remove blank spaces preceding the line count.
      # Remove file name that was returned from "wc".
      # Store line count in "npair" to guide line on velan.

      npair=`wc -l picks.$picknow \
             | sed 's/^  *\(.*\)/\1/' \
             | sed 's/picks.$picknow//' `
      plotline="curve=picks.$picknow npair=$npair curvecolor=white"
    else
      plotline=" "
    fi
    # echo plotline=$plotline

    # plot the cvs
    # cat picks.$picknow
    suximage < cvs.$picknow.su \
	xbox=$XBOXCVS ybox=0 wbox=$WBOXCVS hbox=$HBOX \
        title="CMP $picknow Constant Velocity Stacks" \
        label1=" Time (s)" label2="Velocity (m/s)" \
        f2=$loV d2=$m verbose=0 \
        perc=$myperc n2tic=5 cmap=rgb0  $plotline &

    # if there is a velocity function, display moved out gather
    if [ -s picks.$picknow ]
    then
      # translate picks.$picknow into tnmo/vnmo for sunmo
      sort < picks.$picknow -n |
         mkparfile string1=tnmo string2=vnmo > par.$i
      # apply nmo and plot the moved out gather
      sunmo < panel.$picknow.su par=par.$i cdp=$picknow verbose=0 \
      | $SUXWIGB_OR_XIMAGE  \
         xbox=$XBOXNMOCMP ybox=0 wbox=$WBOXCMP hbox=$HBOX \
         title="CMP $picknow after NMO" \
         label1=" Time (s)" label2="Offset (m)" \
         verbose=0 perc=$myperc key=offset &
    fi

    # compute and plot the semblance/velan
    cat picks.$picknow
    suvelan < panel.$picknow.su nv=$nvs dv=$dvs fv=$fvs \
    | suximage \
	xbox=$XBOXVELAN ybox=0 wbox=$WBOXVELAN hbox=$HBOX perc=99 \
        units="semblance" f2=$fvs d2=$dvs n2tic=5 \
        title="Semblance Plot CMP $picknow" cmap=hsv2 \
        label1=" Time (s)" label2="Velocity (m/s)" \
        legend=1 units=Semblance verbose=0 gridcolor=black \
        grid1=solid grid2=solid \
        mpicks=newpicks.$picknow $plotline

    if [ -s "newpicks.$picknow" ]
    then
      echo "there is a non-zero length newpicks.$picknow file"
      cat newpicks.$picknow
      cp newpicks.$picknow picks.$picknow
    fi

    echo " "
    echo " t-v PICKS CMP $picknow"
    echo "----------------------"
    cat picks.$picknow
    echo "----------------------"
  
    #  rm spanpanel.$picknow.su
    zap xwigb > /dev/null
    zap ximage > /dev/null

  done
  i=`expr $i + 1`

done

#------------------------------------------------
# Create velocity output file
#------------------------------------------------

cdplist=$cmp1

i=2
while [ $i -le $numCMPs ]
do
    eval picknow=\$cmp$i
    cdplist=$cdplist,$picknow
    i=`expr $i + 1`
done
echo cdp=$cdplist \\ >$outpicks

i=1
while [ $i -le $numCMPs ]
do
  sed < par.$i 's/$/ \\/g' >> $outpicks
  i=`expr $i + 1`
done

#------------------------------------------------
# Remove files and exit
#------------------------------------------------
echo " "
echo " The output file of t-v pairs is "$outpicks:
cat $outpicks
rm -f panel.*.su spanpanel.*.su picks.* par.* newpicks.* cvs.*.su
