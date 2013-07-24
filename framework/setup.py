##   Copyright (C) 2009 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Run this from build/framework rather than framework/

import sys, string

try:
    from distutils.core import setup
except ImportError:
    sys.stderr.write('Could not import distutils.\n\n')
    sys.exit(1)

setup(name='madagascar-framework',
      version='1.6-1',
      maintainer='Sergey Fomel',
      maintainer_email='sergey.fomel@gmail.com',
      url='http://www.ahay.org/',
      description='Madagascar Utilities for Reproducible Research',
      py_modules = map(lambda x: 'rsf.'+x,
                       string.split('path doc flow proj prog tex book suproj conf use')),
      scripts=map(lambda x: 'rsf/'+x,
                  string.split('latex2wiki sfdoc sftour')),
      )
