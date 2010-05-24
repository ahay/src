#!/bin/bash

# Makes the current m8r development version into a tarball, removing svn dirs
# Useful for installing on machines that do not have svn.
 
# This is also a prerequisite for making RPMs. Those would be made from a future
# stable release, but for now we need to practice using the development version

SRC=rsfsrc

# Get a fresh copy from the repository, guaranteed not to have local junk in it
svn co https://rsf.svn.sourceforge.net/svnroot/rsf/trunk $SRC

rm -rf $SRC/{admin,user/nobody,api/octave}
find $SRC -name '.svn' -type d -exec rm -rf {} \;

tar -cvzf $SRC.tgz $SRC

rm -rf $SRC
