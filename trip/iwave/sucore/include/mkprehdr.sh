#! /bin/sh
# mkprehdr.sh - make hdr.h from segy.h
# Usage: mkprehdr.sh
# Caution: Don't move the unsigned substitutions below the signed ones!
#
# $Author: john $
# $Source: /usr/local/cwp/src/su/include/RCS/mkprehdr.sh,v $
# $Revision: 1.2 $ ; $Date: 1996/09/09 17:10:28 $

PATH=/bin:/usr/bin

sed '
	/int tracl/,/unass\[/!d
	/;/!d
	s/;.*//
	s/unsigned short/u/
	s/unsigned long/v/
	s/short/h/
	s/long/l/
	s/float/f/
	s/double/z/
	s/int/i/
	s/char/s/
	/unass\[/d
' |
awk '
BEGIN {
printf "static struct {\n\tchar *key;\tchar *type;\tint offs;\n"
printf "} hdr[] = {\n"
}
	{printf "\t{\"%s\",\t\"%s\",\t0},\n", $2, $1}
END {
	printf "};\n"
} '
