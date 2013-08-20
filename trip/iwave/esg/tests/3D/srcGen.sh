#!/bin/sh
rm -f mysrc.su
sx=2300
d=20
fpeak=15
while [ $sx -le 3300 ]
do
  suwaveform type=ricker1 fpeak=$fpeak \
  | sushw key=gelev a=-40 | sushw key=gx a=$sx\
  |sushw key=gy a=500 | sugain scale=100000 >> mysrc.su
  sx=`expr $sx + $d`
done
