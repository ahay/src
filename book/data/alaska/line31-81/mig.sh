#!/bin/bash

cat > temp$$.su

for (( OFFMIN=-5225 ; OFFMIN<5225 ; OFFMIN+=880 ))
do
    OFFMAX=`expr $OFFMIN  + 770`
#    echo $OFFMIN,$OFFMAX
     
    cat temp$$.su | \
    suchw key1=offset key2=sx key3=gx a=0 b=-1 c=1 | \
    suwind key=offset min=$OFFMIN max=$OFFMAX |
    sunmo par=vpick.txt invert=1 | \
    suabshw key=offset | \
    suktmig2d vfile=stkvel1.intp2.bin dx=55
done

rm temp$$.su

