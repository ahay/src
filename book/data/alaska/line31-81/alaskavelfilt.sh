#!/bin/bash
PID=$$

cat >tempin${PID}.su
for (( EP=$1; EP<=$2; EP++))
do
#    echo EP = ${EP} PID = ${PID}
    <tempin${PID}.su suwind \
        key=ep \
        min=${EP} \
        max=${EP} \
    | sustatic hdrs=1 sign=-1 \
    | sugain \
        tpow=1 \
    | sumute \
        xmute=-5280,0,5280 \
	tmute=.488,0,.488 \
    | sugain  agc=1 wagc=.5  \
    > temp${PID}.su
    cat temp${PID}.su | suwind key=offset min=-999998 max=0 \
    | sudipfilt \
        bias=-.003 \
        dx=1 \
        slopes=-.014,-.012,.004,.006 \
        amps=0.0,1.0,1.0,0.0 \
    | suchw \
        key1=offset \
        key2=offset \
        a=0 \
        b=-1


    cat temp${PID}.su | suwind key=offset min=0 max=999998 \
    | sudipfilt \
        bias=.003 \
        dx=1 \
        slopes=-.006,-.004,.012,.014 \
        amps=0.0,1.0,1.0,0.0 

done

rm temp${PID}.su tempin${PID}.su
