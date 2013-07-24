#!/bin/bash
##

# Use: ./vpick2rsf $1 $2 > vpick.rsf
# $1 - where to get picked velocities from, $2 - where to dump column data to

CDPs=(`head -n 1 $1 | sed -e 's/[^0-9]*\([0-9]*\)[^0-9]*/\1 /g'`)

i=0
n2=0
ncdp=${#CDPs[@]}
# CDP sampling
step="55"
# CDP origin
orig="-2612.5"

echo -n > $2

for cdp in `seq 1 $ncdp`; do
    tline=`expr \+ $i \* 2 \+ 2`
    vline=`expr \+ $i \* 2 \+ 3`
    tnmo=`head -n $tline $1 | tail -n 1`
    vnmo=`head -n $vline $1 | tail -n 1`
    Ts=(`echo $tnmo | sed -e 's/[^0-9.]*\([0-9.]*\)[^0-9.]*/\1 /g'`)
    Vs=(`echo $vnmo | sed -e 's/[^0-9.]*\([0-9.]*\)[^0-9.]*/\1 /g'`)
    nt=${#Ts[@]}
    for t in `seq 1 $nt`; do
        echo ${Ts[$t-1]} `echo "${CDPs[$cdp-1]} * $step + $orig" | bc` ${Vs[$t-1]} >> $2
        n2=`expr $n2 \+ 1`
    done
    i=`expr $i \+ 1`
done

echo "data_format=ascii_float"
echo "n1=3"
echo "n2=$n2"
echo "in=$2"
echo

