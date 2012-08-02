#!/bin/sh

awk 'Begin {flag=0;cmp="";done=0} 
{
if ($1 ~ /HANDVEL/ ) {
    if (done==0) {cmp=$2; flag=1 ; done=1 }
    else {print "-1 0 " cmp; cmp=cmp+100; flag=1} 
}else { 
  for( i=1; i<NF; i++) {
     if (flag==1) {print "0.0 " $(i+1)" " cmp; flag=0}
     print $i " " $++i" " cmp
  }
}
}' $1

