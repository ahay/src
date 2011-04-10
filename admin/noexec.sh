#!/bin/bash

# Make non-executable files that should be non-executable

# Copyright (C) 2011 Ioan Vlad
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

find $RSFSRC \
\( \
-name SConstruct -o \
-name SConscript -o \
-name Makefile -o \
-name "*.c" -o \
-name "*.h" -o \
-name "*.f90" -o \
-name "*.f" -o \
-name "*.m" -o \
-name "*.txt" -o \
-name "*.tex" \) \
-perm /u+x,g+x,o+x \
-type f \
-exec svn propdel svn:executable {} \;
