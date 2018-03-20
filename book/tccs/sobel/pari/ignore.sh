#!/bin/bash

ls -1 $1-*.rsf > /dev/null 2>&1
if [ "$?" = "0" ]; then
    temp=$(ls $1-*.rsf)
    for rsf in $temp; do
        data=$( cat $rsf | grep in= | cut -d '"' -f2 | tail -n1)
        rm $rsf
        rm $data
    done
fi
