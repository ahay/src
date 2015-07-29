#!/bin/csh

# Run tests on a subset of "small" examples.

# Copyright (C) 2010 James W. Jennings Jr.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# check $RSFFIGS
if (! $?RSFFIGS) setenv RSFFIGS $RSFROOT/share/madagascar/figs

# list the examples
echo "+++++++ Listing selected examples +++++++"
echo
sfbooklist size=1024 list=filter skipfile=admin/skiplist.txt book

# run the examples
echo "+++++++ Running examples +++++++"
echo
#   need to use "2>&1 >" instead of ">&" because /bin/sh on Ubuntu is /bin/dash
sfbooklist size=1024 list=filter skipfile=admin/skiplist.txt timer=file \
    command="scons > scons.log 2>&1" book

# compare the figures
echo "+++++++ Comparing figures +++++++"
echo
rm book/*/*/*/.rsftest
sfbooklist size=1024 list=filter skipfile=admin/skiplist.txt \
    command="sffiglist rsftest=y list=none > list.log" book

# summarize the tests
echo "+++++++ Summarizing tests +++++++"
echo
sftestlist outfile=testlist.txt book

exit
