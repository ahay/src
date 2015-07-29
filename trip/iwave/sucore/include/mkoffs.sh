#! /bin/sh
# mkoffs.sh - get offset.h for mkhdr.c from segy.h
# Usage: mkoffs.sh
#
# $Author: jkc $
# $Source: /NeXTMount_3.1b/usr/local/cwp/src/su/include/RCS/mkoffs.sh,v $
# $Revision: 1.3 $ ; $Date: 94/07/08 13:48:17 $

PATH=/bin:/usr/bin

cmd=`basename $0`

sed '
	/int tracl/,/unass\[/!d
	/;/!d
	s/;.*//
	/tracl/d
	/unass\[/d
' |
awk '
BEGIN {
	i = 1
}
{
	field = $NF
	printf "\thdr[%d].offs = ", i
	printf "\n\t\toffsetof(segy, %s);\n", field
	++i
}
'
