#! /bin/sh
# File: 010_Velocity_Interp_IVA.sh

#############################################################################
# Velocity interpolation for text file in format for sunmo to prepare it for
# input to suktmig2d
# Credits:
# 2011 Befriko Murdianto - provided script to seisunix@mailman.mines.edu.
# 2011 Schleicher -  Used to process Alaska land line 31-81
#############################################################################


# Input data and directories
VELDIR=.
INVEL=vels.txt
OUTVEL1=stkvel1.intp1.bin
OUTVEL2=stkvel1.intp2.bin

# Define temporal and spatial samples
nsamp=3071
dt=0.004
fsamp=0
ncdpin=10
fcdpin=851
dcdpin=500
ncdpout=4841
fcdpout=838
dcdpout=1

nlines=`grep -v cdp $INVEL | grep -v "#" | wc -l | awk '{print $1}'`
currline=1

>$OUTVEL1

# 1D (fast dimension) velocity interpolation
while [ $currline -le $nlines ]
	do
	nextline=`expr $currline + 1`
	ncolx=`grep -v cdp $INVEL | grep -v "#" | head -$currline | tail -1 | wc -c`
	lcolx=`expr $ncolx - 2`
	ncoly=`grep -v cdp $INVEL | grep -v "#" | head -$nextline | tail -1 | wc -c`
	lcoly=`expr $ncoly - 2`
	Xin=`grep -v cdp $INVEL | grep -v "#" | head -$currline | tail -1 | cut -c6-$lcolx`
	Yin=`grep -v cdp $INVEL | grep -v "#" | head -$nextline | tail -1 | cut -c6-$lcoly`

	unisam xin=$Xin yin=$Yin nout=$nsamp dxout=$dt fxout=$fsamp method=linear >>$OUTVEL1
	currline=`expr $currline + 2`
done
# 2D velocity interpolation
unisam2 nx1=$nsamp dx1=$dt fx1=$fsamp n1=$nsamp d1=$dt f1=$fsamp nx2=$ncdpin dx2=$dcdpin fx2=$fcdpin n2=$ncdpout d2=$dcdpout f2=$fcdpout method1=linear method2=linear <$OUTVEL1 >$OUTVEL2

rm $OUTVEL1

#ximage n1=$nsamp < $OUTVEL2  d1=$dt f2=$fcdpout d2=$dcdpout cmap=rgb11 &


exit

