#!/bin/bash

# Finds and replaces a string with another in files matching pattern
# Pattern and strings are hard-coded below
# Also, be aware that does not preserve permissions. The tmp file gets
# created with the user's default umask

# Copyright (C) 2010 Ioan Vlad
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

for y in `find $RSFSRC -name SConstruct`;
  do 
    sed "s/env\.Docmerge/env\.RSF_Docmerge/g" $y > sed_tmp.asc;
    mv -f sed_tmp.asc $y;
  done

# Matching ".Include" instead of simply "Include" was intended to find 
# occurences of env.Include as well as clones of the environment named 
# differently, i.e. env_f90.Include. A couple of occurences were also found in
# $RSFSRC/configure.py

