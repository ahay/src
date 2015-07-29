#!/bin/bash

# Makes the current m8r development version into a tarball, without svn dirs
# Useful for installing on machines that do not have svn.

# Argument 1: the directory and tarball name. If not given, default is RSFSRC
# Argument 2: SVN version number. If not given, default is current version
# Argument 3: Whether to rm -rf ('clean') or not (no arg given) the downloaded dir


# Argument 1 -- tarball and directory name:
if [ -z $1 ]; then
  SRC='RSFSRC'
else
  SRC=$1
fi

# Argument 2 -- SVN version number:
if [ -z $2 ]; then
  REL='-r '$2
else
  REL=''
fi

# Get a fresh copy from the repository, guaranteed not to have local junk in it
svn export $REL https://rsf.svn.sourceforge.net/svnroot/rsf/trunk $SRC

rm -rf $SRC/{admin,user/nobody,api/octave,scons}

tar -cvzf $SRC.tar.gz $SRC

# Argument 3 -- cleanup:
if test "$3" = "clean"; then
  rm -rf $SRC
fi
